#!/usr/bin/env python3
"""Build cumulative, date-bounded viral protein snapshots for jackhmmer target DBs.

NOTE: the `datasets` CLI has no --released-before flag (only --released-after /
--updated-after). So this pulls the full accession/release-date table for the taxon
ONCE, filters locally per cutoff, and downloads each snapshot by explicit accession list.
"""
import csv, shutil, subprocess, sys, tempfile, time, zipfile
from datetime import datetime
from pathlib import Path

# ---- CONFIG: fill these in before running ----
TAXON = "Orthopoxvirus"                       # e.g. "Orthopoxvirus" or an NCBI Taxonomy ID
CUTOFF_DATES = ["12/31/2015", "12/31/2018", "06/30/2023"]   # e.g. ["12/31/2015", "12/31/2018", "06/30/2023"]
OUTPUT_ROOT = Path("./data/snapshots")
MAX_RETRIES, RETRY_DELAY = 3, 10
# -----------------------------------------------

def log(msg): print(f"[{datetime.now():%Y-%m-%d %H:%M:%S}] {msg}", file=sys.stderr)
def die(msg): log(f"ERROR: {msg}"); sys.exit(1)
def run(cmd, **kw): return subprocess.run(cmd, capture_output=True, text=True, **kw)

def run_with_heartbeat(cmd, label, interval=15):
    """Like run(), but for commands with no native progress output: prints an elapsed-time
    heartbeat every `interval` seconds so a multi-minute call doesn't look hung. Captures to
    temp files rather than pipes since stdout can be several MB (a pipe would deadlock if it
    filled before we read it)."""
    start = time.monotonic()
    with tempfile.TemporaryFile(mode="w+") as out_f, tempfile.TemporaryFile(mode="w+") as err_f:
        proc = subprocess.Popen(cmd, stdout=out_f, stderr=err_f, text=True)
        while proc.poll() is None:
            time.sleep(interval)
            if proc.poll() is None:
                log(f"{label} ({int(time.monotonic() - start)}s elapsed)...")
        out_f.seek(0); stdout = out_f.read()
        err_f.seek(0); stderr = err_f.read()
    return subprocess.CompletedProcess(cmd, proc.returncode, stdout, stderr)

def check_prereqs():
    if not shutil.which("datasets"):
        die("ncbi-datasets-cli not found. Install: conda install -c conda-forge ncbi-datasets-cli")
    if not shutil.which("dataformat"):
        die("dataformat not found (ships with ncbi-datasets-cli).")

def parse_date(s):
    s = s.strip().split("T")[0]
    for fmt in ("%Y-%m-%d", "%m/%d/%Y"):
        try:
            return datetime.strptime(s, fmt)
        except ValueError:
            pass
    return None

def fetch_all_metadata():
    log(f"Pulling full accession/release-date table for taxon '{TAXON}' (one-time)...")
    r1 = run_with_heartbeat(["datasets", "summary", "virus", "genome", "taxon", TAXON, "--as-json-lines"],
                             label="Still pulling accession/release-date table")
    if not r1.stdout.strip():
        die(f"Taxon '{TAXON}' returned no results — check spelling/ID, don't guess.")
    r2 = subprocess.run(["dataformat", "tsv", "virus-genome", "--fields", "accession,release-date"],
                         input=r1.stdout, capture_output=True, text=True)
    lines = [l for l in r2.stdout.strip().splitlines()[1:] if l.strip()]
    records = []
    for line in lines:
        parts = line.split("\t")
        if len(parts) < 2:
            continue
        dt = parse_date(parts[1])
        if dt:
            records.append((parts[0], dt))
    log(f"{len(records)} total records with parseable release dates.")
    return records

def download_by_accessions(acc_file, zip_path):
    # No capture_output here: `datasets download` prints its own live progress bar
    # (bytes downloaded + package validation), which we want visible on screen rather
    # than buffered and thrown away until the process exits.
    delay = RETRY_DELAY
    for attempt in range(1, MAX_RETRIES + 1):
        returncode = subprocess.run(["datasets", "download", "virus", "genome", "accession",
                 "--inputfile", str(acc_file), "--include", "protein",
                 "--filename", str(zip_path)]).returncode
        if returncode == 0:
            return True
        log(f"Download attempt {attempt} failed (see CLI output above for details).")
        time.sleep(delay)
        delay *= 2
    return False

def process_snapshot(cutoff, records, writer):
    tag = cutoff.replace("/", "_")
    d = OUTPUT_ROOT / tag
    d.mkdir(parents=True, exist_ok=True)
    cutoff_dt = parse_date(cutoff)
    if cutoff_dt is None:
        die(f"Could not parse cutoff date '{cutoff}'.")

    matched = [(acc, dt) for acc, dt in records if dt <= cutoff_dt]
    log(f"=== Snapshot {tag}: {len(matched)} accessions released on/before {cutoff} ===")

    if not matched:
        log(f"No accessions on/before {cutoff} — expected for early cutoffs, recording empty snapshot.")
        (d / "protein.faa").touch()
        writer.writerow([tag, cutoff, 0, "", "", str(d), "OK"])
        return

    acc_file = d / "accessions.txt"
    acc_file.write_text("\n".join(acc for acc, _ in matched) + "\n")

    zip_path = d / f"snap_{tag}.zip"
    if not download_by_accessions(acc_file, zip_path):
        log(f"FAILED to download snapshot for {cutoff} after {MAX_RETRIES} attempts.")
        writer.writerow([tag, cutoff, "FAILED", "", "", str(d), "DOWNLOAD_FAILED"])
        return

    with zipfile.ZipFile(zip_path) as z:
        z.extractall(d)

    src, dst = d / "ncbi_dataset" / "data" / "protein.faa", d / "protein.faa"
    n_seq = 0
    if src.exists():
        shutil.move(src, dst)
        n_seq = sum(1 for line in open(dst) if line.startswith(">"))
    else:
        dst.touch()

    zip_path.unlink()
    shutil.rmtree(d / "ncbi_dataset", ignore_errors=True)
    for extra in ("README.md", "md5sum.txt"):
        (d / extra).unlink(missing_ok=True)

    min_dt = min(dt for _, dt in matched)
    max_dt = max(dt for _, dt in matched)
    status = "OK" if max_dt <= cutoff_dt else "QC_FAILED"  # should always hold — we filtered on this ourselves
    if n_seq < len(matched):
        log(f"Note: {len(matched)} accessions requested but only {n_seq} protein sequences returned "
            f"(some records may lack CDS/protein annotation — expected, not necessarily a bug).")
    log(f"Snapshot {tag}: {n_seq} sequences, release-date range [{min_dt.date()} .. {max_dt.date()}]")
    writer.writerow([tag, cutoff, n_seq, min_dt.date().isoformat(), max_dt.date().isoformat(), str(d), status])

def main():
    if TAXON == "<PLACEHOLDER>":
        die("Set TAXON before running.")
    if any("PLACEHOLDER" in c for c in CUTOFF_DATES):
        die("Set CUTOFF_DATES before running.")

    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    check_prereqs()
    records = fetch_all_metadata()

    with open(OUTPUT_ROOT / "summary.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["snapshot_tag", "cutoff_date", "n_sequences", "min_release_date", "max_release_date", "path", "status"])
        for cutoff in CUTOFF_DATES:
            process_snapshot(cutoff, records, writer)

    log(f"Done. Summary at {OUTPUT_ROOT / 'summary.csv'}")

if __name__ == "__main__":
    main()