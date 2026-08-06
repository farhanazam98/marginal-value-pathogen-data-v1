#!/usr/bin/env python3
"""Step 1: submit the query to EBI's jackhmmer web service, poll to completion,
download the alignment of significant hits, and print sanity checks.

API reference (verified live, 2026-07-29, not from memory):
  Submit:   POST https://www.ebi.ac.uk/Tools/hmmer/api/v1/search/jackhmmer
  Poll:     GET  https://www.ebi.ac.uk/Tools/hmmer/api/v1/result/{job_id}
            jackhmmer runs as a *chain* of iteration jobs. While queued this
            endpoint can 500 with an HTML error page instead of returning
            JSON -- that is normal "not ready yet" behavior for this
            endpoint/algo combination, not a fatal error. Once ready it
            returns a JSON array, one entry per iteration completed so far;
            the last entry's own "id" (not the original submission id) is
            the id used for downloads.
  Downloads: POST https://www.ebi.ac.uk/Tools/hmmer/api/v1/download/{iteration_id}/{format}
             triggers async generation (204, no body). Then poll
             GET https://www.ebi.ac.uk/Tools/hmmer/api/v1/download/{iteration_id}
             until the entry for that format has status "AVAILABLE", which
             includes a direct static "url" to the gzipped file.
  Valid database ids (checked live via GET .../search/databases): pfam (hmm),
  refprot, swissprot, uniprot, pdb, rp15, rp35, rp55, rp75, mgnify30_c2,
  mgnify30_c5_fl, mgnify30_c5_ppfam. UniRef100 is NOT offered by this API.
"""

import gzip
import json
import time

import requests
from Bio import SeqIO

QUERY_FASTA = "data/pssm_pipeline/query.fasta"
OUT_ALIGNMENT = "data/pssm_pipeline/msa_raw.sto"
OUT_RAW_RESPONSE = "data/pssm_pipeline/msa_raw_api_response.json"

BASE = "https://www.ebi.ac.uk/Tools/hmmer/api/v1"
DATABASE = "uniprot"  # closest available analog to UniRef100 (see README) -- EBI HMMER API has no UniRef100 option
BITSCORE_PER_RESIDUE = 0.3

POLL_START_DELAY = 10.0
POLL_MAX_DELAY = 60.0
POLL_BACKOFF = 1.3
POLL_TIMEOUT_SECONDS = 60 * 60  # 1 hour ceiling; full UniProt + iteration is slow


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


def main():
    record = next(SeqIO.parse(QUERY_FASTA, "fasta"))
    seq = str(record.seq)
    L = len(seq)
    T = BITSCORE_PER_RESIDUE * L

    print(f"Query length L = {L}")
    print(f"BITSCORE_PER_RESIDUE = {BITSCORE_PER_RESIDUE}")
    print(f"Absolute inclusion threshold T = {BITSCORE_PER_RESIDUE} x {L} = {T} bits")
    print(f"Database = {DATABASE}")

    fasta_text = f">{record.id}\n{seq}\n"
    payload = {
        "input": fasta_text,
        "database": DATABASE,
        "threshold": "bitscore",
        "T": T,
        "incT": T,
        "domT": T,
        "incdomT": T,
    }

    print("\nSubmitting jackhmmer search...")
    submit_resp = requests.post(
        f"{BASE}/search/jackhmmer",
        headers={"Content-Type": "application/json", "Accept": "application/json"},
        json=payload,
    )
    submit_resp.raise_for_status()
    submit_json = submit_resp.json()
    job_id = submit_json["id"]
    print(f"Job ID: {job_id}")

    def is_search_done(parsed):
        if not isinstance(parsed, list) or not parsed:
            return False, "waiting for first iteration"
        last = parsed[-1]
        status = last["status"]
        desc = f"iteration {last['iteration']}, status={status}, convergence_stats={last['convergence_stats']}"
        if status == "SUCCESS":
            return True, desc
        if status in ("FAILURE", "ERROR"):
            raise RuntimeError(f"jackhmmer job failed: {last}")
        return False, desc

    iteration_array = poll_until(lambda: fetch_result(job_id), is_search_done, "search")
    final_iteration = iteration_array[-1]
    iteration_id = final_iteration["id"]
    n_iterations = len(iteration_array)
    print(f"\nConverged/finished after {n_iterations} iteration(s) (server iteration index {final_iteration['iteration']})")
    print(f"Final convergence stats: {final_iteration['convergence_stats']}")

    job_details = requests.get(f"{BASE}/search/{job_id}").json()
    full_result = requests.get(f"{BASE}/result/{iteration_id}", params={"page_size": 10000}).json()
    hits = full_result["result"]["hits"]
    print(f"\nNumber of hits returned: {len(hits)}")

    print("\nTop 10 hit descriptions:")
    for h in hits[:10]:
        desc = h.get("metadata", {}).get("description") if isinstance(h.get("metadata"), dict) else None
        print(f"  score={h['score']:8.1f}  evalue={h['evalue']:.2e}  {h.get('name')}  {desc or ''}")

    print("\nTriggering alignment (Stockholm) generation...")
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

    downloads_list = poll_until(fetch_downloads_list, is_download_ready, "download generation")
    sto_entry = next(d for d in downloads_list if d["format"] == "stockholm")
    print(f"\nDownloading alignment from {sto_entry['url']} ({sto_entry['size']} bytes gzipped)")
    gz_bytes = requests.get(sto_entry["url"]).content
    sto_text = gzip.decompress(gz_bytes).decode("utf-8")
    with open(OUT_ALIGNMENT, "w") as f:
        f.write(sto_text)

    with open(OUT_RAW_RESPONSE, "w") as f:
        json.dump(
            {
                "payload_sent": payload,
                "submit_response": submit_json,
                "job_details": job_details,
                "iteration_array": iteration_array,
                "full_result_stats": full_result["result"]["stats"],
                "downloads_list": downloads_list,
            },
            f,
            indent=2,
        )
    print(f"Wrote raw API response/debug info to {OUT_RAW_RESPONSE}")
    print(f"Wrote alignment to {OUT_ALIGNMENT}")

    print("\nSanity checks:")
    aligned_records = list(SeqIO.parse(OUT_ALIGNMENT, "stockholm"))
    print(f"  Sequences in downloaded alignment: {len(aligned_records)}")

    exact_self_hits = [r for r in aligned_records if str(r.seq).replace("-", "").replace(".", "").upper() == seq.upper()]
    print(f"  Rows in alignment with sequence identical to query: {len(exact_self_hits)}")
    if exact_self_hits:
        print(f"    e.g. {exact_self_hits[0].id}")
    else:
        print("    WARNING: query sequence itself was not found verbatim among the hits.")
        print("    (Top hit below should still be near-identical if this is a known reference protein.)")
        top = hits[0]
        print(f"    Top hit: score={top['score']:.1f} evalue={top['evalue']:.2e} name={top.get('name')}")


if __name__ == "__main__":
    main()
