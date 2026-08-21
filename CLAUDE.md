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

## Pipeline mechanics

- Per-protein settings (protein, bit-score threshold, DMS assays) are read
  from a config file selected via the `PROTEIN_CONFIG` env var (default
  `config/spike.yaml`) — see README's "Configuring which protein" for the
  currently active config.
- Steps 05/06 fan out over every DMS assay listed in the active config,
  producing one row per `(protein, year, assay)` in `data/sweep_results.csv`,
  keyed by `(protein, tag, dms_id)`.
- Sandboxes are keyed by `(protein, year)` at `$SWEEP_ROOT/<protein>/<year>`,
  each with a per-protein PID lock, so different proteins can sweep
  concurrently but a given `(protein, year)` pair cannot run twice at once.
- The bit-score threshold can be swept too, independent of the per-config
  default: `scripts/sweep/run_threshold_sweep.sh` walks a `(year × threshold)`
  grid via a `BITSCORE_PER_RESIDUE` env override (honored in
  `config.load_config()`, so the search and the reuse fingerprint both see
  it), tagging cells `<year>_t<thr>` so they coexist in one
  `sweep_results.csv` — see README's "Running the bit-score threshold sweep".

## Gotchas

- **jackhmmer doesn't parallelize past ~2 cores per job.** Don't increase
  `--cpu` to speed up a single search; get parallelism from running more
  concurrent jobs instead.

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