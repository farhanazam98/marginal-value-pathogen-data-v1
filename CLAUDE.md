# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

See `README.md` for setup, the PSSM pipeline (steps, dependency chain, how to
run them), data acquisition, data layout, and the methodology invariants to
preserve when modifying the pipeline. Keep the two in sync when either
changes — don't duplicate README content here; add to it only what a human
README wouldn't carry (agent-facing notes, in-flight state).

## Current status

**Goal:** get a curve of PSSM's mutation-effect-prediction performance versus
database snapshot year, using whichever UniRef100 snapshots we currently have
on disk. The protein may not stay the current one in the repo (Spike) — we
were originally scoped to Tier A only, but Spike may not have enough related
sequences there, so we may switch to a different protein.

**Notes (as of 2026-08-10):**
- Protein: SARS-CoV-2 Spike, full-length precursor, 1273 aa
  (`data/protein.fasta`, header still generic `>my_protein`) — the one
  candidate to possibly swap out per the goal above.
- Bit-score threshold: still a single hardcoded
  `BITSCORE_PER_RESIDUE = 0.3` in `01_jackhmmer_search.py`, not yet varied
  per protein.
- `run_tier_a_download.sh` ran successfully.
- `run_tier_b_download.sh` overloads the EC2 instance and makes it
  unresponsive — don't just re-run it as-is; the cause needs investigating
  (e.g. concurrency/worker count, resource limits) before Tier B is
  retried.
- Pipeline steps 00-06 run end-to-end (local jackhmmer vs. the local
  UniRef100 2010 snapshot, Spike, threshold 0.3 bits/residue): **Spearman
  rho = 0.175** (95% bootstrap CI [0.142, 0.208], n=3802; excluding the
  1254 imputed variants, rho = 0.203, n=2548). This reproduces the number
  already recorded in `README.md`'s Step 6 to three decimal places, and is
  consistent with the ~0.17 the same local-jackhmmer pipeline produced on
  the Mac checkout — confirming the pipeline is deterministic and gives
  the same result across machines, not machine-dependent. This is the 2010
  point of the curve below.
- **Tier A sweep complete** (2011-2018 run 2026-08-10; 2010 from an earlier
  probe). 9/9 years DONE, no errors. Per-year results in
  `data/sweep_results.csv` (columns documented in
  `data/sweep_results_dictionary.md`), logs in `logs/sweep/`, driver in
  `scripts/sweep/`.
- **rho declines as the snapshot grows**: 0.175 (2010, 4.1 GB) to 0.100
  (2018, 58.8 GB), despite alignment depth rising monotonically
  (Neff 213 to 732). More homologs, worse DMS correlation.
- Caveat: the 95% CIs overlap heavily (2010 [0.142, 0.208] vs 2018
  [0.065, 0.133]) and `imputed_frac` swings 0.17-0.42 across years,
  tracking rho — imputed variants all take one constant value, so they
  enter the Spearman as a tied block. Check `spearman_rho_excl_imputed`
  (a shallower 0.203 to 0.141) before plotting, and state which column any
  chart shows.

This repo is checked out on two machines (a local Mac and an EC2 instance),
each with its own downloaded snapshots and `data/snapshots` layout —
deliberately not committed, since they differ per machine. See
`CLAUDE.local.md` (gitignored, machine-specific — e.g. which database path
`SEQ_DB` should actually point to, and which snapshot years are downloaded
on *this* machine) for the current machine's state.
