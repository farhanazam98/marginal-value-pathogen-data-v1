#!/usr/bin/env python3
"""Fetch UniRef100 FASTA for a set of years from UniProt's previous_major_releases mirror.

Historical releases only ship UniRef100 bundled with UniRef50/90 inside one
combined uniref{YYYY}_{NN}.tar.gz -- no standalone UniRef100 file exists for
any archived year (confirmed against both UniProt and EBI mirrors). This
script streams that archive, extracts only the UniRef100 member, and aborts
the connection the instant that member is fully read, so UniRef50/90 are
never requested.

Archive layout (verified directly against release-2015_01 and release-2023_01
by range-fetching just the tar headers, no full downloads):

    uniref{YYYY}_{NN}.tar.gz            (gzip stream)
     |- tar member #1: uniref100.tar     (always first)
     |    |- README
     |    |- uniref100.release_note
     |    `- uniref100.xml.gz            (the payload)
     |- tar member #2: uniref90.tar
     `- tar member #3: uniref50.tar

The extracted uniref100.xml.gz bytes are piped through a named pipe (FIFO)
into scripts/xml_to_fasta.py, unmodified, running as a subprocess. That
script is already calibrated (14.73 MB/s of compressed input per worker) and
already handles the actual XML->FASTA conversion, so nothing here re-parses
XML. Because it's a FIFO rather than a real file, no compressed or
decompressed XML ever touches disk -- only the final FASTA does.

Concurrency: `workers = 1 if download_MBps <= 6.45 else
ceil(download_MBps / 14.73)`, per the calibration above. Measured once from
the first release's download, then held fixed for the rest of the batch --
network throughput is stable enough within a run that re-measuring
continuously isn't worth the added complexity.

Long-running-download support is deliberately scoped to what this workload
actually needs. Within one release: a dropped connection reconnects with an
HTTP Range request and keeps feeding the same live decompressor, so a network
blip costs a reconnect, not a restart. Across a killed/restarted batch: a
release already on disk with valid stats is skipped, so a restart only
re-pays whatever release was in flight, not the whole multi-hour run. True
mid-file resume across a process restart isn't attempted -- the tar/gzip
decompression state can't be serialized -- but since each release only ever
transfers its own UniRef100 share (about 43% of the combined archive) rather
than the full combined file, that's an acceptable trade rather than something
worth a resumable on-disk staging format.
"""

import argparse
import csv
import json
import math
import os
import subprocess
import sys
import tarfile
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import requests

REPO_ROOT = Path(__file__).resolve().parent.parent
XML_TO_FASTA = Path(__file__).resolve().parent / "xml_to_fasta.py"
ARCHIVE_SIZES_CSV = REPO_ROOT / "data" / "uniprotref_yearly_archive_sizes.csv"

BASE_URL = "https://ftp.uniprot.org/pub/databases/uniprot/previous_major_releases"

# This project's locked-in 13-year set (README, Phase 3: "Year set: Option 2").
OPTION2_YEARS = [2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2020, 2022, 2024, 2026]

# 2010's row in uniprotref_yearly_archive_sizes.csv has no cluster count on
# record (uniref100.release_note wasn't captured for that year); the exact
# figure was verified by hand during xml_to_fasta.py's calibration run
# instead (see its docstring: "9,808,438 vs 9,808,438 published, exact").
CLUSTER_COUNT_FALLBACK = {2010: 9_808_438}

CONNECT_TIMEOUT = 10
READ_TIMEOUT = 30
MAX_RECONNECT_ATTEMPTS = 8


def release_for_year(year):
    """This project's single release-per-year convention: release-{year}_01."""
    return f"{year}_01"


def archive_url(release):
    return f"{BASE_URL}/release-{release}/uniref/uniref{release}.tar.gz"


def parser_workers(download_mbps):
    """Worker count per the calibration in xml_to_fasta.py / the README's Phase 4 rule."""
    if download_mbps <= 6.45:
        return 1
    return math.ceil(download_mbps / 14.73)


class ResilientHTTPStream:
    """A file-like, read()-only view over a URL that reconnects with an HTTP
    Range request (and keeps feeding the same caller-side decompressor) if
    the connection drops mid-transfer, instead of restarting from byte 0."""

    def __init__(self, url, session):
        self.url = url
        self.session = session
        self.pos = 0
        self._response = None
        self._connect()

    def _connect(self):
        headers = {"Range": f"bytes={self.pos}-"} if self.pos else {}
        resp = self.session.get(
            self.url, headers=headers, stream=True,
            timeout=(CONNECT_TIMEOUT, READ_TIMEOUT),
        )
        resp.raise_for_status()
        self._response = resp

    def read(self, size=-1):
        amt = size if size and size > 0 else 1 << 20
        last_err = None
        for attempt in range(MAX_RECONNECT_ATTEMPTS):
            try:
                chunk = self._response.raw.read(amt)
                self.pos += len(chunk)
                return chunk
            except (requests.exceptions.ChunkedEncodingError,
                    requests.exceptions.ConnectionError,
                    requests.exceptions.Timeout,
                    ConnectionError) as exc:
                last_err = exc
                backoff = min(2 ** attempt, 60)
                print(f"  ... connection error at byte {self.pos} ({exc}), "
                      f"retrying in {backoff}s", file=sys.stderr, flush=True)
                time.sleep(backoff)
                self.close()
                self._connect()
        raise RuntimeError(f"exceeded {MAX_RECONNECT_ATTEMPTS} reconnect attempts: {last_err}")

    def close(self):
        if self._response is not None:
            self._response.close()
            self._response = None


def extract_uniref100_xml_gz(stream):
    """Walk the combined archive to the uniref100.xml.gz member, streaming.

    Returns a sequential file-like object positioned at the start of that
    member's (still gzip-compressed) bytes. Reading it to EOF consumes
    exactly that member -- nothing from uniref90/50 is ever pulled off the
    wire, since those tar members come later in the same ordered stream and
    are simply never read.
    """
    outer = tarfile.open(fileobj=stream, mode="r|gz")
    first = outer.next()
    if first is None or first.name != "uniref100.tar":
        raise RuntimeError(
            f"expected uniref100.tar as the first archive member, got {first and first.name!r}"
        )
    inner_stream = outer.extractfile(first)

    inner = tarfile.open(fileobj=inner_stream, mode="r|")
    for member in inner:
        if member.name == "uniref100.xml.gz":
            return inner.extractfile(member)
    raise RuntimeError("uniref100.xml.gz not found inside uniref100.tar")


def load_cluster_counts():
    counts = dict(CLUSTER_COUNT_FALLBACK)
    with open(ARCHIVE_SIZES_CSV, newline="") as f:
        for row in csv.DictReader(f):
            value = row["uniref100_cluster_count"]
            if value and value != "NA":
                counts[int(row["year"])] = int(value)
    return counts


def integrity_ok(year, stats_path, cluster_counts):
    if not stats_path.exists():
        return False
    try:
        stats = json.loads(stats_path.read_text())
    except (OSError, json.JSONDecodeError):
        return False
    expected = cluster_counts.get(year)
    if expected is None:
        return bool(stats.get("n_seqs"))
    return stats.get("n_seqs") == expected


def download_release(year, output_dir, cluster_counts):
    """Download, extract, and convert one year's UniRef100. Returns (year, download_MBps)."""
    release = release_for_year(year)
    url = archive_url(release)
    fasta_path = output_dir / f"uniref100_{release}.fasta"
    stats_path = output_dir / f"uniref100_{release}.stats.json"

    if integrity_ok(year, stats_path, cluster_counts):
        print(f"[{year}] already complete, skipping ({fasta_path})", flush=True)
        return year, None

    print(f"[{year}] fetching {url}", flush=True)
    t0 = time.time()
    session = requests.Session()
    stream = ResilientHTTPStream(url, session)

    with tempfile.TemporaryDirectory(prefix=f"uniref100_{release}_") as tmp_dir:
        fifo_path = Path(tmp_dir) / "uniref100.xml.gz"
        os.mkfifo(fifo_path)

        proc = subprocess.Popen([
            sys.executable, str(XML_TO_FASTA),
            str(fifo_path), str(fasta_path),
            "--stats-json", str(stats_path),
        ])

        try:
            xml_gz_stream = extract_uniref100_xml_gz(stream)
            with open(fifo_path, "wb") as fifo:
                while True:
                    chunk = xml_gz_stream.read(1 << 20)
                    if not chunk:
                        break
                    fifo.write(chunk)
        finally:
            # Early abort: we stop reading the moment uniref100.xml.gz is
            # fully consumed, so uniref90/50 -- later in the same ordered
            # stream -- are never requested.
            stream.close()

        returncode = proc.wait()

    elapsed = time.time() - t0
    download_mbps = (stream.pos / (1 << 20)) / elapsed if elapsed > 0 else 0.0

    if returncode != 0 or not integrity_ok(year, stats_path, cluster_counts):
        print(f"[{year}] FAILED (xml_to_fasta.py exit {returncode}); "
              f"removing partial output", file=sys.stderr, flush=True)
        fasta_path.unlink(missing_ok=True)
        stats_path.unlink(missing_ok=True)
        raise RuntimeError(f"release {release} failed integrity/conversion check")

    print(f"[{year}] done in {elapsed:.0f}s, "
          f"{stream.pos / (1 << 20):.0f} MB fetched at {download_mbps:.2f} MB/s "
          f"-> {fasta_path}", flush=True)
    return year, download_mbps


def run_batch(years, output_dir, max_workers_override):
    output_dir.mkdir(parents=True, exist_ok=True)
    cluster_counts = load_cluster_counts()

    remaining = list(years)
    failures = []

    if max_workers_override is None and len(remaining) > 1:
        # Calibrate worker count off the first release's real throughput,
        # then hold it fixed for the rest of the batch.
        first_year = remaining.pop(0)
        try:
            _, mbps = download_release(first_year, output_dir, cluster_counts)
        except RuntimeError as exc:
            print(f"[{first_year}] {exc}", file=sys.stderr, flush=True)
            failures.append(first_year)
            mbps = None
        worker_count = parser_workers(mbps) if mbps else 1
        print(f"calibrated {worker_count} worker(s) from {first_year} "
              f"({mbps:.2f} MB/s)" if mbps else
              f"calibration unavailable, defaulting to 1 worker", flush=True)
    else:
        worker_count = max_workers_override or 1

    if remaining:
        with ThreadPoolExecutor(max_workers=worker_count) as pool:
            futures = {
                pool.submit(download_release, year, output_dir, cluster_counts): year
                for year in remaining
            }
            for future in as_completed(futures):
                year = futures[future]
                try:
                    future.result()
                except RuntimeError as exc:
                    print(f"[{year}] {exc}", file=sys.stderr, flush=True)
                    failures.append(year)

    if failures:
        print(f"\n{len(failures)} release(s) failed: {sorted(failures)}", file=sys.stderr)
        sys.exit(1)


def parse_years(raw_years):
    if raw_years == ["option2"]:
        return OPTION2_YEARS
    return sorted(int(y) for y in raw_years)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--years", nargs="+", required=True,
                     help="years to fetch, e.g. --years 2010 2011 2015, "
                          "or --years option2 for this project's locked-in 13-year set")
    ap.add_argument("--output-dir", default=str(REPO_ROOT / "data"),
                     help="where uniref100_{release}.fasta and .stats.json are written")
    ap.add_argument("--max-workers", type=int, default=None,
                     help="override the calibrated worker/concurrency count")
    args = ap.parse_args()

    years = parse_years(args.years)
    run_batch(years, Path(args.output_dir), args.max_workers)


if __name__ == "__main__":
    main()
