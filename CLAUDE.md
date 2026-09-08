# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

See `README.md` for setup, the PSSM pipeline (steps, dependency chain, how to
run them), data acquisition, data layout, current status (active config,
sweep progress, findings), and the methodology invariants to preserve when
modifying the pipeline. Keep the two in sync when either changes — don't
duplicate README content here; add to it only what a human README wouldn't
carry (agent-facing notes, in-flight state).

## Working conventions

**Keep changes as simple as reasonably possible.** Prefer the smallest change
that does the job. Don't add abstraction, configuration, options, or new files
that aren't needed yet. When editing docs, prefer cutting over reorganizing,
and don't restate the same fact in two places — put it where a reader would
look for it and link to it from elsewhere.

**Docs describe current state, not history.** Don't reference old or removed
behavior ("no longer a constant", "used to use X") — describe what the code
does now. Git history and GitHub issues carry the before; don't keep a
running changelog or strikethrough TODO list in this file or the README.

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

**Match EVEREST's methodology when changing the pipeline.** This project
reimplements the alignment-based half of EVEREST (Gurev/Youssef/Marks,
bioRxiv 2025.08.04.668549, <https://doi.org/10.1101/2025.08.04.668549>). When a
design choice isn't
already pinned down by README's "Key methodology to preserve" section, default
to whatever EVEREST does rather than inventing a new approach, and add it to
that section if you deviate.

## Pipeline mechanics

- Per-protein settings (protein, bit-score threshold, DMS assays) are read
  from a config file selected via the `PROTEIN_CONFIG` env var (default
  `config/spike.yaml`) — see README's "Configuring which protein" for the
  currently active config.
- Steps 05/06 fan out over every DMS assay listed in the active config,
  producing one row per `(protein, year, assay)` in `data/sweep_results.csv`,
  keyed by `(protein, tag, dms_id)`.
- Sandboxes are keyed by `(protein, cell)` at `$SWEEP_ROOT/<protein>/<cell>`,
  each with a per-protein PID lock, so different proteins can sweep
  concurrently but a given cell cannot run twice at once. Every cell is named
  `<year>_t<thr>`, so a cell's threshold is readable from its name.
- The bit-score threshold can be swept too, independent of the per-config
  default: `scripts/sweep/run_threshold_sweep.sh` walks a `(year × threshold)`
  grid via a `BITSCORE_PER_RESIDUE` env override (honored in
  `config.load_config()`, so the search and the reuse fingerprint both see
  it), tagging cells `<year>_t<thr>` so they coexist in one
  `sweep_results.csv` — see README's "Running the bit-score threshold sweep".
  - With no `-t`, the driver sweeps EVEREST's exact grid
    {0.01, 0.03, 0.05, 0.1, 0.3, 0.5} by default (the low end is the regime a
    sparse family needs to clear the depth floor). This makes an ordinary run 6
    sequential DB scans (~6× a single year sweep, since the PID lock serializes
    thresholds); pass `-t` to narrow it. In practice the lowest thresholds can
    be dropped per protein when those scans exceed the wall-clock budget on the
    larger snapshots, so a protein's swept grid is often a subset. The grids
    actually committed are: spike {0.1, 0.2, 0.3, 0.4, 0.5}; flu_h1_ha
    {0.05, 0.1, 0.3, 0.5}; hiv_env, dengue_polg and protease
    {0.03, 0.04, 0.05, 0.1, 0.3, 0.5}.
  - Alignment selection (Neff/L > 1, then max fraction ≥90% ID) lives only in
    `plot_threshold_sweep.py`; the scoring pipeline emits the inputs (`Neff`,
    `Neff_at_90pct_identity`) but no step picks an alignment.

## Gotchas

- **Never run `collect.py` or commit `data/sweep_results.csv` on a per-protein
  sweep branch.** That CSV is the one artifact not partitioned by protein, so a
  regenerated copy on the branch makes two proteins' branches collide on merge.
  Commit only `data/sweep/<protein>` metas; rebuild the CSV once at combine, on
  the full-snapshot machine — see README's "Running proteins across separate
  machines".
- **jackhmmer doesn't parallelize past ~2 cores per job.** Don't increase
  `--cpu` to speed up a single search; get parallelism from running more
  concurrent jobs instead.
- **Bit-score threshold changes are free.** In calibration, a 5x threshold
  change (0.1 to 0.5 bits/residue) moved wall-clock time by only 2.5% — the
  full database scan happens regardless of threshold, so tuning it per
  protein, or running the threshold sweep, costs nothing extra. Memory is
  similarly a non-issue at every stage (peak RSS stayed under 150 MB
  regardless of database size).
- **Storage throughput, not CPU, is the EC2 risk.** Search is I/O-bound at
  ~385 MB/s per job in calibration, above the default gp3 EBS baseline
  (125 MB/s) — a year's FASTA needs to be page-cached or on local NVMe, or
  search runs roughly 3x slower than measured. See `docs/calibration.csv` for the
  raw numbers behind these three.

## Known issues and decisions

- **`prop_90` selects on the wrong quantity — not being fixed.** It should mean
  "fraction of the alignment within 90% identity *of the query*" (the query row
  of the identity matrix) but computes how self-redundant the alignment is with
  itself, which favours small, collapsed, near-duplicate alignments. Sites:
  `03_weights.py:70-96`, `plot_threshold_sweep.py:119` (the ratio) and `:73`
  (`idxmax` over it).
  **Decision:** don't repair the selector. Report every swept threshold, and mark
  the best-performing cell per `(protein, year)`. That cell is an *oracle* — it
  is picked using the DMS being predicted, so it is a theoretical upper bound on
  what threshold choice could buy, not a score the pipeline can reach on its own.
  Always label it as a ceiling; it is defensible as a bound and indefensible as a
  result. Any "more data doesn't help" framing waits on this all-threshold
  analysis.

- **Every tracked symlink under `data/sweep` is an absolute EC2 path** — 1622 of
  them, all `/home/ec2-user/...`, so they dangle on any clone. Sandbox
  scaffolding, not results: `collect.py` reads only the JSON metas and `STATUS`.
  Fix: stop tracking them, or emit them relative in `run_year.sh`.

- **Open decision: how to make the year axis consistent.** Spike ran 13 years;
  the other four ran 14, having added a 2025 snapshot as an EVEREST comparison
  point. So no figure can currently put all five on one axis without a caveat.
  The options, none picked yet: drop the 2025 cells and standardise on the
  shared 13 years (2010–2018, 2020, 2022, 2024, 2026); sweep spike's 2025 to
  bring it up to 14; or keep the asymmetry and state which axis every figure
  uses. Decide before any cross-protein figure ships.

## Verifying changes

After any change to the scoring pipeline (steps 04–06), re-run the sandbox
for a single already-completed year (e.g. 2018) and confirm `spearman_rho`
in `data/sweep_results.csv` matches the last-committed value for that
`(protein, year, assay)` row to 3 decimal places. A changed value means the
change altered scoring behavior, not just style — flag it before continuing
rather than assuming it's an improvement.

This repo is checked out on two machines (a local Mac and an EC2 instance),
each with its own downloaded snapshots and `data/snapshots` layout —
deliberately not committed, since they differ per machine. See
`CLAUDE.local.md` (gitignored, machine-specific — e.g. which database path
`SEQ_DB` should actually point to, and which snapshot years are downloaded
on *this* machine) for the current machine's state.
