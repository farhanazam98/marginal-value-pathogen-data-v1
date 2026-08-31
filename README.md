# marginal-value-pathogen-data-v1

## Project

A research project measuring whether PSSM mutation-effect prediction accuracy
(Spearman rho vs. a DMS assay) changes as UniRef100 database snapshots grow
2010→2026, accounting for sequence diversity (e.g. Neff@90%ID). It is a minimal
reimplementation of the alignment-based half of the EVEREST pipeline
(Gurev/Youssef/Marks, bioRxiv 2025.08.04.668549; local copy `docs/EVEREST.pdf`).
Two proteins are swept so far: SARS-CoV-2 Spike (vs. Starr 2020 DMS) and
SARS-CoV-2 main protease (vs. Flynn fitness DMS) — see Current status.

## Results

![PSSM accuracy vs. UniRef100 snapshot year, 2010-2026](plots/pssm_accuracy_vs_snapshot_year.png)
![Protease PSSM accuracy vs. UniRef100 snapshot year, 2010-2026](plots/protease_accuracy_vs_snapshot_year.png)

Spike's accuracy (Spearman's rho) drops as the database grows through 2018, then
holds roughly flat through 2026; protease's does not decline at all — see
Current status → Findings for the numbers and caveats. Source data:
`data/sweep_results.csv` (columns documented in
`data/sweep_results_dictionary.md`); regenerate a plot with
`PROTEIN_CONFIG=config/<name>.yaml python scripts/sweep/plot.py`.

## Setup

```bash
./setup.sh
conda activate marginal-value-pathogen-data
```

`setup.sh` installs Miniconda if missing, then creates/updates the
`marginal-value-pathogen-data` conda env from `environment.yml` (Python 3.11 plus
matplotlib, numpy, pandas, scipy, requests, biopython, and HMMER 3.4, which
provides the `jackhmmer` binary). Re-running it is safe. Pipeline steps that need
scipy (e.g. `06_evaluate.py`) require this env — the system `python3` won't do.

No test suite, linter, or CI: correctness is checked by inline per-step sanity
checks, described in each pipeline script and the sections below.

## Running the PSSM pipeline

Steps live under `scripts/pssm_pipeline/`, `00_*.py`–`06_*.py`, one per stage,
each runnable standalone from the repo root and each writing its checkpoint(s)
to `data/pssm_pipeline/` (gitignored, regenerable):

```bash
python scripts/pssm_pipeline/00_setup.py             # canonicalize + validate query.fasta
python scripts/pssm_pipeline/01_jackhmmer_search.py  # jackhmmer search -> raw MSA
python scripts/pssm_pipeline/02_clean_msa.py         # Stockholm -> N x L matrix, EVEREST A.3.1 filters
python scripts/pssm_pipeline/03_weights.py           # 99%-identity sequence reweighting, Neff
python scripts/pssm_pipeline/04_pssm.py              # per-column aa frequency table
python scripts/pssm_pipeline/05_score.py             # log-odds score every DMS variant, per assay
python scripts/pssm_pipeline/06_evaluate.py          # Spearman rho vs each DMS assay, bootstrap CI, scatter
```

Each step consumes only the previous step's checkpoint, so any one can be re-run
after a change: `00 query.fasta` → `01 msa_raw.sto` → `02 msa_clean.npy` →
`03 weights.npy` → `04 pssm.npy` → `05 predictions.csv` → `06 rho + scatter.png`.

Step 01 needs a local UniRef100 FASTA snapshot on disk (see Data acquisition).
Point it at one with the `SEQ_DB` env var (default `data/uniref100_2010.fasta`),
e.g. `SEQ_DB=data/snapshots/uniref100_2015_01.fasta`. The bit-score threshold
comes from the active config's `bitscore_per_residue` (see Configuring which
protein); it is length-normalized (bits/residue × query length), not e-value
based, so hit stringency stays constant as snapshots grow.

## Configuring which protein

Which protein the pipeline runs is set by a per-protein YAML config, selected
with `PROTEIN_CONFIG` (default `config/spike.yaml`):

```yaml
name: SARS2_Spike
query_fasta: data/proteins/protein.fasta
bitscore_per_residue: 0.3
assays:
  - {id: starr_binding,    csv: data/dms/SARS2_RBD_Starr_binding_dms.csv, label: Starr 2020 ACE2 binding}
  - {id: starr_expression, csv: data/dms/SARS2_RBD_Starr_expression.csv,  label: Starr 2020 RBD expression}
```

Steps 00–01 read the query and threshold; steps 05–06 score *every* assay in
`assays` against the single PSSM, writing `predictions_<id>.csv`,
`scatter_<id>.png`, and per-assay metas. Building the MSA/PSSM (00–04) is the
expensive per-protein work; adding an assay on top is a cheap fan-out.

To add a protein, copy `config/spike.yaml`, point it at that protein's query
FASTA and DMS CSV(s), and set `PROTEIN_CONFIG`. Each DMS must use the same
residue numbering as its query — step 05 checks the WT residue in every `mutant`
against the query and stops on any mismatch; there is no coordinate-offset
support.

## Running the sweep across years

`scripts/sweep/run_sweep.sh` runs steps 00–06 once per snapshot year:

```bash
PROTEIN_CONFIG=config/spike.yaml scripts/sweep/run_sweep.sh -j 6 2010 2011 2012 2013 2014 2015 2016 2017 2018
python scripts/sweep/collect.py
```

Each year builds one PSSM and scores every assay against it, in its own sandbox
at `$SWEEP_ROOT/<protein>/<year>` (`SWEEP_ROOT` default `data/sweep/`,
git-tracked): symlinks to the shared scripts, config, and inputs, plus its own
`data/pssm_pipeline/`. `<protein>` is the config filename stem and the PID lock
is per-protein, so two proteins can sweep at once — keep the combined `-j` within
the machine's cores. `-j` caps concurrent years (default 6); jackhmmer only uses
~2 cores per job, so higher just oversubscribes.

Reruns are cheap: a year whose sandbox already holds a PSSM for the same query +
threshold skips the search (00–04) and only re-scores (05–06). Changing the query
or threshold forces a rebuild.

`collect.py` rebuilds `data/sweep_results.csv` (one row per `(protein, year,
assay)`, columns in `data/sweep_results_dictionary.md`) from whatever sandboxes
exist under `$SWEEP_ROOT` — safe to run mid-sweep. `plot.py` draws one line per
assay. **Run `collect.py` / `plot.py` only on the machine with the full snapshot
set** — on a machine missing years they write a partial or mixed-schema table
over the good one. The sweep itself only touches sandbox checkpoints, never these
two tracked files. Undo an accidental regen with `git checkout --
data/sweep_results.csv plots/pssm_accuracy_vs_snapshot_year.png`.

## Running the bit-score threshold sweep

`scripts/sweep/run_threshold_sweep.sh` runs a 2D `(year × threshold)` grid
instead of holding the threshold at the config value — to check whether a
stricter or looser threshold changes the rho-vs-year shape:

```bash
scripts/sweep/run_threshold_sweep.sh -j 6 2010 2011 ... 2026
python scripts/sweep/collect.py
```

It is a thin wrapper over `run_sweep.sh`: for each threshold it sets the
`BITSCORE_PER_RESIDUE` env override (read by `config.load_config()`, so both the
search and the PSSM-reuse fingerprint see it) and runs a year sweep whose cells
are tagged `<year>_t<thr>` (e.g. `2018_t0.05`), each in its own sandbox
`data/sweep/<protein>/<year>_t<thr>`, all merging into one `sweep_results.csv`
with the threshold in the `bitscore_per_residue` / `threshold_bits` columns.

With no `-t` it sweeps EVEREST's grid `{0.01 0.03 0.05 0.1 0.3 0.5}`; pass
`-t "..."` to narrow it. Thresholds run **sequentially** (the per-protein PID
lock forbids concurrent sweeps under one root), years concurrently within a
threshold, so the full grid is ~6× a single year sweep — six independent
database scans, not one reused. The default grid includes the spike baseline
0.3, whose `_t0.3` cells re-derive and (the pipeline being deterministic)
validate the plain year sweep; keep 0.3 in any `-t` to preserve that check. Do
**not** set `SWEEP_ROOT` — every cell must land in the default `data/sweep` root
or `collect.py` silently truncates the merged table (its only tell is the
`Wrote … (N rows)` count).

`plot.py` stays the single-threshold rho-vs-year curve (config baseline threshold
only); the threshold-vs-rho comparison is a separate figure.

## Running proteins across separate machines

To sweep several proteins at once, give each its own machine and git branch, then
combine at the end. Everything a run writes is already keyed by protein
(sandboxes, PID lock, `logs/`, the per-protein plot) — except
`data/sweep_results.csv`, so the one rule is: **don't regenerate that file on a
per-protein machine; rebuild it once at the end.**

Per protein, on its own machine (an image/clone already carrying the repo and
this protein's snapshots):

```bash
git checkout main && git pull            # pick up newly-committed protein configs
git checkout -b sweep/<protein>
PROTEIN_CONFIG=config/<protein>.yaml scripts/sweep/run_threshold_sweep.sh -j 6 \
  2010 2011 2012 2013 2014 2015 2016 2017 2018 2020 2022 2024 2026
git add data/sweep/<protein>             # commit only this protein's metas + STATUS
git commit -m "Sweep <protein> across the threshold grid"
git push -u origin sweep/<protein>
```

Do **not** run `collect.py` or `plot.py` here — a rewritten
`data/sweep_results.csv` reaching the branch would collide when two proteins'
branches merge. Check `git status` doesn't list it before pushing.

Combine once, on the full-snapshot machine (so `db_n_seqs` / `db_n_residues` fill
from its `.stats.json` sidecars), after merging every branch:

```bash
git checkout main
git merge sweep/<A> sweep/<B> sweep/<C>  # disjoint meta dirs, no conflict
python scripts/sweep/collect.py          # rebuild the full CSV from every protein's metas
for p in <A> <B> <C>; do PROTEIN_CONFIG=config/$p.yaml python scripts/sweep/plot.py; done
git add data/sweep_results.csv plots/ && git commit -m "Combine protein sweeps"
```

A new protein needs its `config/<protein>.yaml`, query FASTA, and DMS CSV
committed before launch. On an image-cloned instance the multi-GB snapshots ride
along on the cloned volume — confirm they are present and non-zero first.

## Data acquisition (UniRef100 snapshots)

`scripts/download_uniref100.py` fetches one year's UniRef100 FASTA at a time:

```bash
python scripts/download_uniref100.py --years option2   # this project's locked 13-year set
python scripts/download_uniref100.py --years 2010 2015 --output-dir data/snapshots
```

Archived UniProt releases only ship a combined `uniref{YYYY}_{NN}.tar.gz`
(UniRef50+90+100 together). The script streams that tarball, extracts only the
`uniref100.xml.gz` member (always first in the stream), and drops the connection
right after.

Two mirrors carry these at an identical path, `ftp.uniprot.org` and
`ftp.ebi.ac.uk`, and either can silently go slow or hang at the TLS handshake. At
the start of a batch the script reads 8 MB from each and pins the *faster* one
(not the first to answer — a mirror can send headers promptly then deliver at
0.4 MB/s) for every year in the batch, since resuming a dropped transfer with an
HTTP `Range` request into a live decompressor can't switch hosts mid-file. A
mirror dying mid-batch fails those years; rerun to re-probe and resume from disk.

`--min-free-gb` (default 150) refuses to start a year when the output volume is
below that floor; peak disk exceeds the final FASTA because a year's `.xml.gz` is
kept until its FASTA passes the integrity check. Download and parse are separate
concurrent phases: `--download-workers` (default 4, network/disk-bound) and
`--parse-workers` (default 3, CPU-bound, capped low for predictable memory).
`scripts/xml_to_fasta.py` is a standalone streaming parser (no DOM).

`scripts/run_tier_a_download.sh` and `scripts/run_tier_b_download.sh` launch the
2010–2018 and 2020/2022/2024/2026 batches as detached `nohup` jobs (PID file
under `data/snapshots/`, logs in `logs/`) that survive an SSH disconnect; each
prints its `kill -0` / `kill` commands on start. `data/snapshots` is expected to
be a symlink to scratch/NVMe — the root volume has no room for these. Which years
are on disk differs per machine and isn't committed; check `ls data/snapshots`.

## Data layout

- `config/*.yaml` — per-protein configs (query FASTA, threshold, assay list);
  `config/spike.yaml` is the default, selected by `PROTEIN_CONFIG`. Tracked.
- `data/proteins/` — one query FASTA per protein (`protein.fasta` = full-length
  Spike, 1273 aa; `protease_protein.fasta` = Mpro, 306 aa). Tracked.
- `data/dms/` — one CSV per assay: `SARS2_RBD_Starr_binding_dms.csv` and
  `SARS2_RBD_Starr_expression.csv` (Starr 2020 Spike, RBD positions 331–531);
  `SARS2_MRPO_Flynn_dms.csv` (protease fitness). Each shares its query's
  numbering. Tracked.
- `data/pssm_pipeline/` — gitignored checkpoint dir for a by-hand pipeline run
  from the repo root; steps 05–06 write one set per assay.
- `data/sweep/<protein>/<cell>/` — per-cell sweep sandboxes (`<cell>` is a year,
  or `<year>_t<thr>`), each a self-contained pipeline working dir. Only the JSON
  metas + `STATUS` are tracked — enough for `collect.py` to rebuild
  `sweep_results.csv`; the heavy binaries are gitignored.
- `data/snapshots/` — gitignored multi-GB UniRef100 FASTA snapshots plus
  `.stats.json` sidecars; usually a symlink to scratch/NVMe.
- `data/uniprotref_yearly_archive_sizes.csv` — combined UniRef50+90+100 archive
  size per year; feeds the growth plot and `download_uniref100.py`'s expected
  cluster counts.
- `data/sweep_results.csv` — one row per `(protein, year, assay)`, keyed by
  `(protein, tag, dms_id)`, carrying every step's metrics; columns in
  `data/sweep_results_dictionary.md`. Produced by `collect.py`, not hand-edited.
- `plots/` — tracked PNGs linked from this README, regenerated by `plot.py` (one
  per `PROTEIN_CONFIG`). `expression_vs_binding.png` is a standalone figure, not
  produced by any current script.
- `calibration.csv` — timed measurements from early calibration runs; the
  findings drawn from it live in CLAUDE.md's Gotchas, not read by any script.

## Current status

**Goal:** a curve of PSSM mutation-effect-prediction accuracy versus database
snapshot year, for whichever protein is under study.

### Configured proteins

Per-protein settings are read from the `PROTEIN_CONFIG` config (see Configuring
which protein). Two proteins are fully swept:

- **Spike** (`config/spike.yaml`) — SARS-CoV-2 Spike, 1273 aa. Threshold `0.3`
  bits/residue. Assays: `starr_binding`, `starr_expression` (both Starr 2020).
- **Protease / Mpro** (`config/protease.yaml`) — SARS-CoV-2 Mpro, 306 aa.
  Threshold `0.1` bits/residue. Assay: `flynn_fitness`.

Spike's `0.3` was picked by maximizing rho on Spike itself — validation leakage.
Protease runs at a fixed `0.1` baseline (one of EVEREST's six swept thresholds),
not rho-tuned. Spike has since run the full `(year × threshold)` grid and
threshold turns out not to matter for its finding (see Findings), but protease
has only run at `0.1`, so the Spike-vs-protease contrast is still confounded by
threshold as well as protein identity. The final round re-sweeps protease across
the same grid to close this.

### Final-round proteins (configured, sweep pending)

Four more proteins, picked to span how much each virus's sequence data has grown
since 2010 — so the rho-vs-year shape can be read against a slow-evolving control
and a hypervariable case, not just one virus. Each has a config, query FASTA, and
EVEREST DMS CSV(s) committed; none is swept. Run each on its own
snapshot-bearing machine through the full threshold grid.

| Virus | Protein | Config | DMS (EVEREST) | aa | Role |
|-------|---------|--------|---------------|----|------|
| Influenza A H1N1 | Hemagglutinin (H1 HA) | `flu_h1_ha` | `IAV_H1_HA_Doud`, `IAV_H1_HA_Wu` | 565 | influenza; two labs share one PSSM |
| HIV-1 (BG505) | Envelope (Env) | `hiv_env` | `HIV1_BG505_ENV_Haddox` | 860 | HIV |
| AAV2 | Capsid (VP1) | `aav2_capsid` | `AAV2_CAPSD_Sinai` | 735 | control virus |
| Dengue | Polyprotein region (POLG) | `dengue_polg` | `DENV_POLG_Suphatrakul` | 900 | priority arbovirus |

Each config's `bitscore_per_residue` is `0.1` as a baseline, but the final round
runs the full `(year × threshold)` grid via `run_threshold_sweep.sh`, which
overrides it.

### Progress

- **Spike sweep complete**: 13 years (2010–2018, 2020, 2022, 2024, 2026), both
  assays — 26 rows in `data/sweep_results.csv`.
- **Protease sweep complete**: the same 13 years, one assay — 13 rows.
- Sandboxes live under `data/sweep/<protein>/<year>/`; only the JSON metas +
  `STATUS` are tracked, enough to rebuild the CSV but not to regenerate
  alignments/PSSMs without rerunning against the snapshots.

### Findings

- **Spike: rho declines, then plateaus, as the snapshot grows.** For
  `starr_binding`, rho falls from 0.175 (2010, 4.1 GB) to 0.100 (2018, 58.8 GB),
  then holds at ~0.10–0.12 across 2020–2026 (up to 219 GB) — more homologs, never
  better agreement with the DMS data. `starr_expression` traces the same shape at
  a higher baseline (0.248 → 0.172, then ~0.17–0.20).
  - Endpoint 95% bootstrap CIs for `starr_binding` are disjoint (2010
    [0.142, 0.208] vs 2018 [0.065, 0.133]), so the 2010→2018 decline is real;
    2016 breaks the trend upward, and the post-2018 years all overlap, so the
    plateau is genuinely flat rather than a smooth curve.
  - Alignment depth `Neff_over_L` (homologs per column, corrected for
    near-duplicates) climbs monotonically (0.25 → 1.89) and only crosses
    EVEREST's depth-adequacy floor of 1.0 in 2020 — every 2010–2018 year would
    fail that check, which limits how much the absolute rho values in that range
    can bear.
  - `imputed_frac` (share of DMS variants whose alignment column got dropped, so
    they get a constant fill instead of a real prediction) swings 0.16–0.42,
    tracking `L_final` (816–879 of 1273 columns surviving).
  - `jackhmmer_converged` is `False` for 5 of the 13 years (2010, 2011, 2013,
    2024, 2026) — the 5-round cap was hit while still finding ~one new hit per
    round, not a failure.

- **Protease: rho is much higher, and does not decline.** `flynn_fitness` rho
  starts at 0.541 (2010) and drifts up to 0.572 (2024) / 0.567 (2026) — the
  opposite direction from Spike. The drift is shallow and the CIs mostly overlap,
  so read it as "no decline," not a confirmed increase.
  - `jackhmmer` converges cleanly in every one of the 13 years, unlike Spike.
  - `imputed_frac` stays near zero (0–0.013): `L_final` is 302–306 of the
    protein's 306 columns every year.
  - `Neff_over_L` only crosses the depth floor in 2022, yet the earlier
    depth-inadequate years still give rho on par with the rest — unlike Spike,
    where crossing the floor lines up with the point rho stopped declining.
  - Protease differs from Spike in both bit-score threshold and protein identity
    (shorter, more conserved), so this contrast is suggestive rather than
    controlled — it doesn't say whether the threshold or the protein explains the
    difference.

- **Bit-score threshold is not load-bearing for the Spike finding.** Across
  {0.1–0.5} bits/residue the decline-then-plateau shape holds; only 0.5 raises
  rho, and only by returning degenerate near-duplicate-only alignments that fail
  the depth floor. See `plots/spike_threshold_sweep*.png`.

> Update this section as status changes; keep CLAUDE.md pointing here rather than
> duplicating it. Completed work lives in git history, not a TODO list here.

## Key methodology to preserve when modifying the pipeline

These choices follow the EVEREST paper (`docs/EVEREST.pdf`); where a decision
isn't pinned down below, default to whatever EVEREST does and note any deviation
here.

- **Bit-score thresholds, not e-values**, for the jackhmmer search: an e-value
  cutoff would silently tighten as later-year snapshots grow, while bits/residue
  × query length keeps hit stringency constant across the 2010–2026 sweep.
- **Sequence reweighting at 99% identity** (theta=0.01), not the more common 80%,
  because two sequences differing by 1% can already have meaningfully different
  fitness for this protein. Neff/L (effective sequences per column) is the
  depth-adequacy check against EVEREST's floor of 1.0.
- **Alignment selection is DMS-blind and per protein** (EVEREST Methods A.6.1):
  among alignments with Neff/L > 1, use the one with the highest *fraction* of
  sequences within 90% identity of the query (`Neff@90%ID / Neff`). Never select
  on rho — using the DMS you're predicting to choose the pipeline is validation
  leakage. Sweep the threshold per protein and **keep scoring every swept
  threshold, not just the selected one**: the selection rule can pick a
  degenerate alignment (on Spike, 2026 selects 0.5 — 80% imputed, rho 0.03), and
  only the full set surfaces that.
- **Column/sequence filtering in `02_clean_msa.py`** (drop columns >50% gaps,
  drop sequences <50% query coverage) is computed against the *original* query
  positions independently for both filters, not against each other's output —
  preserve this independence if the filters are ever touched.
- **Coordinate mapping is offset-free by construction** (single-query jackhmmer
  profile ⇒ match column *i* is query position *i*; each DMS and the query
  already share numbering) — step 05 reconciles this per assay and stops if it
  fails to hold. If the query or search strategy changes such that this breaks,
  every downstream step's coordinate assumptions need re-verification.
