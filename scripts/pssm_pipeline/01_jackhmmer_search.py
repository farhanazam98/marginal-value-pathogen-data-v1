#!/usr/bin/env python3
"""Step 1 (revised): sweep jackhmmer over EVEREST's inclusion thresholds.

Instead of one search at a single threshold, build one raw candidate alignment
per EVEREST length-normalized bit-score threshold (Methods A.4.1). Every search
uses the same EBI jackhmmer web service against the same `uniprot` database that
produced the rho=0.248 benchmark -- only the inclusion threshold varies, so any
change downstream is attributable to threshold alone (apples-to-apples). The API
offers no UniRef100 (verified live, printed below), so UniRef100 is a separate,
later, local-jackhmmer track, not part of this sweep.

The 0.3 threshold was already run; its output is reused by renaming, not
re-queried. The remaining five thresholds are searched here, serially, by
looping the same submit -> poll -> download flow used before. Each candidate is
checkpointed to its own labeled file, so any one can be inspected without a
rerun and an interrupted sweep resumes without repeating completed searches.

API reference (verified live, 2026-07-29, not from memory):
  Submit:   POST https://www.ebi.ac.uk/Tools/hmmer/api/v1/search/jackhmmer
  Poll:     GET  https://www.ebi.ac.uk/Tools/hmmer/api/v1/result/{job_id}
            jackhmmer runs as a *chain* of iteration jobs. While queued this
            endpoint can 500 with an HTML error page instead of returning
            JSON -- that is normal "not ready yet" behavior, not fatal. Once
            ready it returns a JSON array, one entry per iteration completed so
            far; the last entry's own "id" (not the submission id) is the id
            used for downloads.
  Downloads: POST https://www.ebi.ac.uk/Tools/hmmer/api/v1/download/{iteration_id}/{format}
             triggers async generation (204, no body). Then poll
             GET https://www.ebi.ac.uk/Tools/hmmer/api/v1/download/{iteration_id}
             until the entry for that format has status "AVAILABLE", which
             includes a direct static "url" to the gzipped file.
  Databases (checked live via GET .../search/databases): pfam, refprot,
  swissprot, uniprot, pdb, rp15, rp35, rp55, rp75, mgnify30_c2, mgnify30_c5_fl,
  mgnify30_c5_ppfam. UniRef100 is NOT offered by this API.
"""

import gzip
import json
import os
import time

import requests
from Bio import SeqIO

QUERY_FASTA = "data/pssm_pipeline/query.fasta"
OUT_DIR = "data/pssm_pipeline"

BASE = "https://www.ebi.ac.uk/Tools/hmmer/api/v1"
DATABASE = "uniprot"  # closest available analog to UniRef100 (see README) -- EBI HMMER API has no UniRef100 option

# EVEREST Methods A.4.1 length-normalized bit-score thresholds. 0.3 is already
# on disk (reused via rename); the rest are searched by this script.
REUSED_THRESHOLD = 0.3
THRESHOLDS_TO_SEARCH = [0.5, 0.1, 0.05, 0.03, 0.01]
ALL_THRESHOLDS = sorted([REUSED_THRESHOLD] + THRESHOLDS_TO_SEARCH, reverse=True)

LEGACY_RAW = f"{OUT_DIR}/msa_raw.sto"  # the pre-sweep single-search output (threshold 0.3)
LEGACY_RESPONSE = f"{OUT_DIR}/msa_raw_api_response.json"

POLL_START_DELAY = 10.0
POLL_MAX_DELAY = 60.0
POLL_BACKOFF = 1.3
POLL_TIMEOUT_SECONDS = 60 * 60  # 1 hour ceiling per search; full UniProt + iteration is slow


def raw_path(threshold):
    return f"{OUT_DIR}/msa_raw_t{threshold}.sto"


def response_path(threshold):
    return f"{OUT_DIR}/msa_raw_api_response_t{threshold}.json"


def poll_until(fetch_fn, is_done_fn, label):
    """Poll fetch_fn() with backoff until is_done_fn(response) is True. Returns the final parsed response."""
    start = time.monotonic()
    delay = POLL_START_DELAY
    while True:
        elapsed = time.monotonic() - start
        if elapsed > POLL_TIMEOUT_SECONDS:
            raise TimeoutError(f"{label}: exceeded {POLL_TIMEOUT_SECONDS}s without completing")
        parsed = fetch_fn()
        if parsed is not None:
            done, status_desc = is_done_fn(parsed)
            print(f"  [{elapsed:6.0f}s] {label}: {status_desc}")
            if done:
                return parsed
        else:
            print(f"  [{elapsed:6.0f}s] {label}: not ready yet")
        time.sleep(delay)
        delay = min(delay * POLL_BACKOFF, POLL_MAX_DELAY)


def fetch_result(job_id):
    resp = requests.get(f"{BASE}/result/{job_id}", headers={"Accept": "application/json"})
    if resp.status_code != 200:
        return None
    try:
        return resp.json()
    except ValueError:
        return None


def print_database_list():
    """Print the live jackhmmer database list and assert UniRef100 is absent."""
    resp = requests.get(f"{BASE}/search/databases", headers={"Accept": "application/json"})
    resp.raise_for_status()
    dbs = resp.json()
    ids = [d["id"] for d in dbs]
    print("Live jackhmmer databases (GET /search/databases):")
    print(f"  {', '.join(ids)}")
    has_uniref = any("uniref" in i.lower() for i in ids)
    print(f"  UniRef100 offered by this API: {has_uniref} "
          f"(sweep uses '{DATABASE}'; UniRef100 is a later local-jackhmmer track)")
    if DATABASE not in ids:
        raise RuntimeError(f"Target database '{DATABASE}' is not in the live list {ids} -- stopping.")


def run_search(threshold, seq, record_id):
    """Submit one jackhmmer search at `threshold`, poll to completion, download the
    Stockholm alignment to raw_path(threshold), and return per-candidate stats."""
    L = len(seq)
    T = threshold * L
    print(f"\n=== threshold {threshold} bits/residue  (absolute T = {threshold} x {L} = {T} bits, db={DATABASE}) ===")

    payload = {
        "input": f">{record_id}\n{seq}\n",
        "database": DATABASE,
        "threshold": "bitscore",
        "T": T,
        "incT": T,
        "domT": T,
        "incdomT": T,
    }

    submit_resp = requests.post(
        f"{BASE}/search/jackhmmer",
        headers={"Content-Type": "application/json", "Accept": "application/json"},
        json=payload,
    )
    submit_resp.raise_for_status()
    submit_json = submit_resp.json()
    job_id = submit_json["id"]
    print(f"  Job ID: {job_id}")

    def is_search_done(parsed):
        if not isinstance(parsed, list) or not parsed:
            return False, "waiting for first iteration"
        last = parsed[-1]
        status = last["status"]
        desc = f"iteration {last['iteration']}, status={status}, convergence_stats={last['convergence_stats']}"
        if status == "SUCCESS":
            return True, desc
        if status in ("FAILURE", "ERROR"):
            raise RuntimeError(f"jackhmmer job failed at threshold {threshold}: {last}")
        return False, desc

    iteration_array = poll_until(lambda: fetch_result(job_id), is_search_done, f"search t{threshold}")
    final_iteration = iteration_array[-1]
    iteration_id = final_iteration["id"]
    n_iterations = len(iteration_array)

    full_result = requests.get(f"{BASE}/result/{iteration_id}", params={"page_size": 10000}).json()
    hits = full_result["result"]["hits"]
    print(f"  Converged after {n_iterations} iteration(s); {len(hits)} hits")

    gen_resp = requests.post(f"{BASE}/download/{iteration_id}/stockholm")
    if gen_resp.status_code not in (200, 204):
        print(f"  WARNING: generate request returned {gen_resp.status_code}: {gen_resp.text[:200]}")

    def fetch_downloads_list():
        resp = requests.get(f"{BASE}/download/{iteration_id}")
        if resp.status_code != 200:
            return None
        return resp.json()

    def is_download_ready(parsed):
        entry = next((d for d in parsed if d["format"] == "stockholm"), None)
        if entry is None:
            return False, "stockholm format entry not present yet"
        return entry["status"] == "AVAILABLE", f"stockholm status={entry['status']}"

    downloads_list = poll_until(fetch_downloads_list, is_download_ready, f"download t{threshold}")
    sto_entry = next(d for d in downloads_list if d["format"] == "stockholm")
    print(f"  Downloading alignment ({sto_entry['size']} bytes gzipped) -> {raw_path(threshold)}")
    gz_bytes = requests.get(sto_entry["url"]).content
    sto_text = gzip.decompress(gz_bytes).decode("utf-8")
    with open(raw_path(threshold), "w") as f:
        f.write(sto_text)

    with open(response_path(threshold), "w") as f:
        json.dump(
            {
                "threshold": threshold,
                "payload_sent": payload,
                "submit_response": submit_json,
                "iteration_array": iteration_array,
                "full_result_stats": full_result["result"]["stats"],
                "downloads_list": downloads_list,
            },
            f,
            indent=2,
        )

    return {"threshold": threshold, "n_iterations": n_iterations, "n_hits": len(hits), "status": "searched"}


def ensure_t03_reused():
    """Reuse the existing 0.3 search by renaming its output to the labeled name.
    Do not re-query. Loud-stop if the 0.3 candidate is genuinely absent."""
    target = raw_path(REUSED_THRESHOLD)
    if os.path.exists(target):
        print(f"Reused 0.3 candidate already present: {target}")
        return
    if os.path.exists(LEGACY_RAW):
        os.rename(LEGACY_RAW, target)
        print(f"Renamed existing 0.3 alignment: {LEGACY_RAW} -> {target} (not re-queried)")
        if os.path.exists(LEGACY_RESPONSE):
            os.rename(LEGACY_RESPONSE, response_path(REUSED_THRESHOLD))
            print(f"Renamed its API response:      {LEGACY_RESPONSE} -> {response_path(REUSED_THRESHOLD)}")
        return
    raise RuntimeError(
        f"Neither {target} nor {LEGACY_RAW} exists -- the reused 0.3 candidate is missing. "
        "Stopping rather than silently searching it fresh (would break the reuse assumption). "
        "Add 0.3 to THRESHOLDS_TO_SEARCH deliberately if you want it re-queried."
    )


def alignment_row_count(path):
    """Count sequence rows in a Stockholm alignment (comparable 'size' per candidate)."""
    return sum(1 for _ in SeqIO.parse(path, "stockholm"))


def main():
    record = next(SeqIO.parse(QUERY_FASTA, "fasta"))
    seq = str(record.seq)
    print(f"Query {record.id}, length L = {len(seq)}\n")

    print_database_list()

    print("\n--- Reuse the already-completed 0.3 search (rename, no re-query) ---")
    ensure_t03_reused()

    print("\n--- Search the remaining thresholds (serial; skips any already on disk) ---")
    searched_stats = {}
    for threshold in THRESHOLDS_TO_SEARCH:
        if os.path.exists(raw_path(threshold)):
            print(f"\n=== threshold {threshold}: {raw_path(threshold)} already present, skipping search ===")
            searched_stats[threshold] = {"threshold": threshold, "status": "already_on_disk"}
            continue
        searched_stats[threshold] = run_search(threshold, seq, record.id)

    # ---------------- Sanity table: one row per candidate ----------------
    print("\n=== Step 1 sanity: candidates present ===")
    print(f"{'threshold':>10}  {'database':>9}  {'status':>14}  {'align_rows':>10}  {'iterations':>10}")
    all_present = True
    for threshold in ALL_THRESHOLDS:
        path = raw_path(threshold)
        if not os.path.exists(path):
            all_present = False
            print(f"{threshold:>10}  {DATABASE:>9}  {'MISSING':>14}  {'-':>10}  {'-':>10}")
            continue
        rows = alignment_row_count(path)
        stat = searched_stats.get(threshold, {})
        if threshold == REUSED_THRESHOLD:
            status = "reused_0.3"
            iters = "n/a (reused)"
        else:
            status = stat.get("status", "?")
            iters = stat.get("n_iterations", "n/a" if status == "already_on_disk" else "?")
        print(f"{threshold:>10}  {DATABASE:>9}  {status:>14}  {rows:>10}  {str(iters):>10}")

    print()
    if not all_present:
        raise RuntimeError("At least one candidate alignment is missing -- see MISSING rows above. Stopping.")
    print(f"All {len(ALL_THRESHOLDS)} candidate alignments present under {OUT_DIR}/msa_raw_t*.sto.")
    print("Ready for Step 2 (clean each candidate).")


if __name__ == "__main__":
    main()
