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
- `run_tier_b_download.sh` downloads and parses each year with independently-
  sized worker pools — see README's Data acquisition section.
- **Tier A sweep complete** (2011-2018 run 2026-08-10; 2010 from an earlier
  probe). 9/9 years DONE, no errors. Per-year results in
  `data/sweep_results.csv`, logs in `logs/sweep/`, driver in `scripts/sweep/`
  (see README's "Running the sweep across years").
- **Tier B download complete** (2026-08-13, all 4 years: 2020/2022/2024/2026,
  ~621 GB). Zero length mismatches, same integrity check as Tier A. All 13
  years of the locked set are now on disk.
- **Tier B sweep complete** (2026-08-14, all 4 years: 2020/2022/2024/2026).
  4/4 DONE.

**Findings:**
- **rho declines then plateaus as the snapshot grows**: 0.175 (2010,
  4.1 GB) falls to 0.100 (2018, 58.8 GB), then holds flat at ~0.10-0.12
  across all four Tier B years (2020-2026, up to 219 GB), while alignment
  depth climbs monotonically the whole way (Neff 213 to 1664). More
  homologs, never better DMS correlation: marginal value is negative
  2010→2018, then zero.
  - Endpoint CIs are disjoint (2010 [0.142, 0.208] vs 2018 [0.065, 0.133]),
    so the 2010→2018 decline is real, but adjacent years overlap and 2016
    breaks the trend upward — a trend, not a smooth curve. The Tier B
    points all sit inside each other's CIs, so the plateau is genuinely
    flat, not a second trend.
  - 2024 and 2026 did not converge within jackhmmer's 5 rounds
    (`jackhmmer_converged=False`).
  - `imputed_frac` swings 0.16-0.42 across years. Imputed variants all
    take one constant value, so they enter the Spearman as a tied block
    carrying no rank information. The fraction is set by `L_final`
    (872-879 in the low-imputation years, 816-844 in the high ones).

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
4. ~~Run the sweep against the 4 Tier B years (2020/2022/2024/2026).~~ Done
   2026-08-14 — 4/4 DONE, full 2010-2026 curve now in
   `data/sweep_results.csv`.
5. Diagnose the rho decline: run
  `scripts/diagnostics/rbd_gap_diagnostic.py` against each year's
  `msa_raw.sto` to check whether later, larger snapshots pull in more
  divergent homologs that gap out specifically in the RBD (positions
  361-413) — which would explain the 50%-gap-column filter dropping more
  of the RBD as snapshots grow, one candidate mechanism for the decline
  alongside the imputed_frac caveat above.
6. Add the alignment-threshold heuristic from the EVEREST paper.

This repo is checked out on two machines (a local Mac and an EC2 instance),
each with its own downloaded snapshots and `data/snapshots` layout —
deliberately not committed, since they differ per machine. See
`CLAUDE.local.md` (gitignored, machine-specific — e.g. which database path
`SEQ_DB` should actually point to, and which snapshot years are downloaded
on *this* machine) for the current machine's state.