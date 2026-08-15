# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

See `README.md` for setup, the PSSM pipeline (steps, dependency chain, how to
run them), data acquisition, data layout, and the methodology invariants to
preserve when modifying the pipeline. Keep the two in sync when either
changes — don't duplicate README content here; add to it only what a human
README wouldn't carry (agent-facing notes, in-flight state).

## Working conventions

**Keep changes as simple as reasonably possible.** Prefer the smallest change
that does the job. Don't add abstraction, configuration, options, or new files
that aren't needed yet. When editing docs, prefer cutting over reorganizing,
and don't restate the same fact in two places — put it where a reader would
look for it and link to it from elsewhere.

**Explain esoteric concepts in plain terms.** Assume the reader has no strong
biology background. Domain jargon (PSSM, MSA, `Neff`, bit-score threshold,
DMS, imputation, homolog) gets a plain-English gloss at first use, and prose
should say what a number *means*, not just report it. Prefer a concrete
example over a definition where one fits.

**Verify the user's understanding before committing anything, especially
prose.** Do not commit until they have confirmed they follow what changed and
why. Walk through the change in plain terms and wait for a response — an
absent objection is not confirmation. This applies to documentation and
written explanations as much as to code.

## Current status

**Goal:** get a curve of PSSM's mutation-effect-prediction performance versus
database snapshot year, using whichever UniRef100 snapshots we currently have
on disk. The protein may not stay the current one in the repo (Spike).

**Configuration (as of 2026-08-12):**
- Protein: SARS-CoV-2 Spike, full-length precursor, 1273 aa
  (`data/protein.fasta`, header still generic `>my_protein`) — the one
  candidate to possibly swap out per the goal above.
- Bit-score threshold: still a single hardcoded
  `BITSCORE_PER_RESIDUE = 0.3` in `01_jackhmmer_search.py`, not yet varied
  per protein.

**Progress:**
- `run_tier_b_download.sh` overloaded the EC2 instance's RAM and made it
  unresponsive. Root cause: `download_uniref100.py` used to fuse downloading
  and parsing into one worker per year (FIFO + subprocess), with the worker
  count picked to match parser throughput — so one knob controlled both
  network and parsing concurrency, and each concurrent worker held an HTTP
  stream, two live levels of tar/gzip decompression, and a parser process at
  once, all for the entire multi-GB transfer. Fixed by separating download
  (writes `.xml.gz` to disk) from parse (reads that file, writes FASTA,
  deletes the `.xml.gz`) into two independently-sized worker pools:
  `--download-workers` (default 4) and `--parse-workers` (default 3, capped
  low since it's the CPU/memory-bound side). Not yet re-tested against an
  actual Tier B run — see Open TODOs.
- **Full 13-year sweep complete** (Tier A 2010-2018 run 2026-08-10,
  reproduced against redownloaded snapshots 2026-08-13; Tier B
  2020/2022/2024/2026 completed 2026-08-14). 13/13 years DONE, no errors.
  Per-year results in `data/sweep_results.csv` (columns documented in
  `data/sweep_results_dictionary.md`; regenerate with
  `scripts/sweep/collect.py`), logs in `logs/sweep/`, driver in
  `scripts/sweep/`. Plotted in `pssm_accuracy_vs_snapshot_year.png`
  (`scripts/sweep/plot.py`) — see README's Results section for the chart
  itself; this section covers the numbers and caveats behind it.

**Findings:**
- **rho by year**: 0.175 (2010, 4.1 GB) → 0.0996 (2018, 58.8 GB) → 0.1230
  (2020, 94.3 GB) → 0.1063 (2022, 132.8 GB) → 0.1067 (2024, 175.8 GB) →
  0.1098 (2026, 218.6 GB). Neff: 213 (2010) → 732 (2018) → 1664 (2026).
  - CIs: 2010 [0.142, 0.208], 2018 [0.065, 0.133], 2026 [0.078, 0.141].
    2010 vs. 2018 disjoint; 2018 vs. 2026 overlap.
  - `imputed_frac` (fraction of DMS variants landing on a column the
    50%-gap filter dropped, so they get a constant fill value instead of a
    prediction): 0.33 (2010), 0.41 (2018), 0.41 (2020), 0.165 (2022), 0.165
    (2024), 0.160 (2026). Full 13-point range: 0.16-0.42, same range as the
    original 9-point set.

**Open TODOs:**
1. ~~Separate the parsing component from the download script.~~ Done — see
   Progress above.
2. ~~Repoint `data/snapshots` at the new EBS volume and redownload Tier A
   there.~~ Done — symlinks to `/data/snapshots`, all 9 years (2010-2018,
   ~188 GB) redownloaded 2026-08-12/13, stats show no length mismatches.
   See `CLAUDE.local.md`.
3. ~~Rerun the pipeline after the previous step and reproduce results from
   before.~~ Done 2026-08-13 — all 9 years (2010-2018) re-run against the
   redownloaded snapshots; every value in `data/sweep_results.csv` matches
   the 2026-08-10 run exactly (including `spearman_rho` for all nine years)
   except wall-clock timing, confirming the pipeline is deterministic and
   the rho-decline finding isn't an artifact of one run.
4. Diagnose the rho decline: run
  `scripts/diagnostics/rbd_gap_diagnostic.py` against each year's
  `msa_raw.sto` to check whether later, larger snapshots pull in more
  divergent homologs that gap out specifically in the RBD (positions
  361-413) — which would explain the 50%-gap-column filter dropping more
  of the RBD as snapshots grow, one candidate mechanism for the decline
  alongside the imputed_frac caveat above.
5. Add the alignment-threshold heuristic from the EVEREST paper.

This repo is checked out on two machines (a local Mac and an EC2 instance),
each with its own downloaded snapshots and `data/snapshots` layout —
deliberately not committed, since they differ per machine. See
`CLAUDE.local.md` (gitignored, machine-specific — e.g. which database path
`SEQ_DB` should actually point to, and which snapshot years are downloaded
on *this* machine) for the current machine's state.