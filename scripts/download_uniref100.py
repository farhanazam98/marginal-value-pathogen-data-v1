#!/usr/bin/env python3
"""Fetch UniRef100 XML for a set of years from UniProt's previous_major_releases mirror.

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

Download and parse are two separate, independently-sized phases connected by
the extracted uniref100.xml.gz file on disk (not a live pipe): downloading a
year writes that file, then a bounded pool of parse workers converts it to
FASTA via scripts/xml_to_fasta.py and deletes the .xml.gz once the FASTA is
verified. --download-workers and --parse-workers size the two phases
independently (see README.md's Data acquisition section for why they used to
be fused, and why that was a problem).

Long-running-download support is deliberately scoped to what this workload
actually needs. Within one release: a dropped connection reconnects with an
HTTP Range request and keeps feeding the same live decompressor, so a network
blip costs a reconnect, not a restart. Across a killed/restarted batch: a
release already fully parsed is skipped, and a release already downloaded but
not yet parsed skips straight to parsing -- so a restart only re-pays whatever
was actually in flight.
"""

import argparse
import csv
import json
import subprocess
import sys
import tarfile
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import requests
import urllib3

REPO_ROOT = Path(__file__).resolve().parent.parent
XML_TO_FASTA = Path(__file__).resolve().parent / "xml_to_fasta.py"
ARCHIVE_SIZES_CSV = REPO_ROOT / "data" / "uniprotref_yearly_archive_sizes.csv"

# ftp.uniprot.org itself hangs at the TLS handshake from this instance
# (confirmed 2026-08-12: TCP connects, ClientHello sent, no ServerHello --
# times out); EBI mirrors the same archives at an identical path layout and
# is reachable.
BASE_URL = "https://ftp.ebi.ac.uk/pub/databases/uniprot/previous_major_releases"

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

# Reading via response.raw (rather than requests' iter_content) means urllib3
# exceptions reach us unwrapped -- requests only translates urllib3 errors
# into requests.exceptions.* at the iter_content layer. A bare
# urllib3.exceptions.ProtocolError (e.g. "Connection reset by peer" mid-tar)
# is what actually killed a prior run: it isn't a subclass of the builtin
# ConnectionError, so it slipped past a narrower catch tuple. HTTPError is
# the common base for ProtocolError/IncompleteRead/ReadTimeoutError/etc.
TRANSIENT_NETWORK_ERRORS = (
    requests.exceptions.ChunkedEncodingError,
    requests.exceptions.ConnectionError,
    requests.exceptions.Timeout,
    urllib3.exceptions.HTTPError,
    ConnectionError,
)


def release_for_year(year):
    """This project's single release-per-year convention: release-{year}_01."""
    return f"{year}_01"


def archive_url(release):
    return f"{BASE_URL}/release-{release}/uniref/uniref{release}.tar.gz"


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
            except TRANSIENT_NETWORK_ERRORS as exc:
                last_err = exc
                backoff = min(2 ** attempt, 60)
                print(f"  ... connection error at byte {self.pos} ({exc}), "
                      f"retrying in {backoff}s", file=sys.stderr, flush=True)
                time.sleep(backoff)
                self.close()
                try:
                    self._connect()
                except TRANSIENT_NETWORK_ERRORS as reconnect_exc:
                    # The reconnect itself hit the same kind of transient
                    # error -- count it as this attempt and let the loop
                    # retry with backoff, rather than propagating and
                    # killing the whole batch.
                    last_err = reconnect_exc
                    continue
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


def download_release(year, output_dir):
    """Download one year's UniRef100 XML (still gzip-compressed) to disk.

    Streams to a .partial path and renames atomically on completion, so a
    half-written file left behind by a crash is never mistaken for a
    complete one.
    """
    release = release_for_year(year)
    url = archive_url(release)
    xml_gz_path = output_dir / f"uniref100_{release}.xml.gz"
    partial_path = xml_gz_path.with_suffix(xml_gz_path.suffix + ".partial")

    print(f"[{year}] fetching {url}", flush=True)
    t0 = time.time()
    session = requests.Session()
    stream = ResilientHTTPStream(url, session)

    try:
        xml_gz_stream = extract_uniref100_xml_gz(stream)
        with open(partial_path, "wb") as f:
            while True:
                chunk = xml_gz_stream.read(1 << 20)
                if not chunk:
                    break
                f.write(chunk)
    except Exception:
        partial_path.unlink(missing_ok=True)
        raise
    finally:
        # Early abort: we stop reading the moment uniref100.xml.gz is fully
        # consumed, so uniref90/50 -- later in the same ordered stream -- are
        # never requested.
        stream.close()

    partial_path.rename(xml_gz_path)

    elapsed = time.time() - t0
    mbps = (stream.pos / (1 << 20)) / elapsed if elapsed > 0 else 0.0
    print(f"[{year}] downloaded in {elapsed:.0f}s, "
          f"{stream.pos / (1 << 20):.0f} MB at {mbps:.2f} MB/s -> {xml_gz_path}",
          flush=True)


def parse_release(year, xml_gz_path, output_dir, cluster_counts):
    """Convert one year's downloaded uniref100.xml.gz to FASTA via
    xml_to_fasta.py, verify it against the expected cluster count, and
    delete the .xml.gz once the FASTA is confirmed good."""
    release = release_for_year(year)
    fasta_path = output_dir / f"uniref100_{release}.fasta"
    stats_path = output_dir / f"uniref100_{release}.stats.json"

    print(f"[{year}] parsing {xml_gz_path}", flush=True)
    result = subprocess.run([
        sys.executable, str(XML_TO_FASTA),
        str(xml_gz_path), str(fasta_path),
        "--stats-json", str(stats_path),
    ])

    if result.returncode != 0 or not integrity_ok(year, stats_path, cluster_counts):
        print(f"[{year}] FAILED (xml_to_fasta.py exit {result.returncode}); "
              f"removing partial output", file=sys.stderr, flush=True)
        fasta_path.unlink(missing_ok=True)
        stats_path.unlink(missing_ok=True)
        raise RuntimeError(f"release {release} failed parsing/integrity check")

    xml_gz_path.unlink()
    print(f"[{year}] parsed -> {fasta_path}", flush=True)


def process_year(year, output_dir, cluster_counts, parse_pool):
    """Get one year to a verified FASTA: skip if already done, download if
    the .xml.gz isn't on disk yet, then parse.

    Runs as a download-pool task, but hands the actual parsing off to the
    separately-sized parse_pool and blocks on its result -- so this thread
    sits idle waiting rather than parsing itself, and parse concurrency
    stays capped at --parse-workers no matter how many downloads are
    in flight at once.
    """
    release = release_for_year(year)
    fasta_path = output_dir / f"uniref100_{release}.fasta"
    stats_path = output_dir / f"uniref100_{release}.stats.json"
    xml_gz_path = output_dir / f"uniref100_{release}.xml.gz"

    if integrity_ok(year, stats_path, cluster_counts):
        print(f"[{year}] already parsed, skipping ({fasta_path})", flush=True)
        return

    if not xml_gz_path.exists():
        download_release(year, output_dir)

    parse_pool.submit(parse_release, year, xml_gz_path, output_dir, cluster_counts).result()


def run_batch(years, output_dir, download_workers, parse_workers):
    output_dir.mkdir(parents=True, exist_ok=True)
    cluster_counts = load_cluster_counts()
    failures = []

    with ThreadPoolExecutor(max_workers=download_workers) as dl_pool, \
         ThreadPoolExecutor(max_workers=parse_workers) as parse_pool:
        futures = {
            dl_pool.submit(process_year, year, output_dir, cluster_counts, parse_pool): year
            for year in years
        }
        for fut in as_completed(futures):
            year = futures[fut]
            try:
                fut.result()
            except Exception as exc:
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
    ap.add_argument("--download-workers", type=int, default=4,
                     help="concurrent downloads -- network/disk-bound, cheap per worker (default 4)")
    ap.add_argument("--parse-workers", type=int, default=3,
                     help="concurrent xml_to_fasta.py conversions -- CPU-bound, capped to bound "
                          "total memory (default 3)")
    args = ap.parse_args()

    years = parse_years(args.years)
    run_batch(years, Path(args.output_dir), args.download_workers, args.parse_workers)


if __name__ == "__main__":
    main()
