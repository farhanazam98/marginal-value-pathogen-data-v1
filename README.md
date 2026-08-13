# marginal-value-pathogen-data-v1

## Project

A research project measuring whether PSSM mutation-effect prediction accuracy
(Spearman rho vs. Starr 2020 DMS) change as UniRef100 database snapshots grow
2010→2026, accounting for sequence diversity (e.g Neff@90%ID). It is a
minimal reimplementation of the alignment-based half of the EVEREST pipeline
(Gurev/Youssef/Marks, bioRxiv 2025.08.04.668549).

## Setup and environment

```bash
./setup.sh
conda activate marginal-value-pathogen-data
```

`setup.sh` installs Miniconda if missing, then creates/updates the
`marginal-value-pathogen-data` conda env from `environment.yml` (Python 3.11,
matplotlib, numpy, pandas, scipy, requests, biopython, and HMMER 3.4 — the
`jackhmmer` binary used by the alignment pipeline). Re-running it is safe.
Some pipeline steps (e.g. `06_evaluate.py`, which needs scipy) require this
conda env specifically — the system `python3` is not sufficient.

There is no test suite, linter, or CI config in this repo. Correctness is
validated per-step via inline sanity checks (see each pipeline script and its
README section) rather than automated tests.

## Running the PSSM pipeline

Steps live under `scripts/pssm_pipeline/`, numbered `00_*.py`–`06_*.py`, one
per stage, each runnable standalone from the repo root and each writing its
checkpoint(s) to `data/pssm_pipeline/` (gitignored — regenerable from the
scripts) so any step's output can be inspected without rerunning earlier ones:

```bash
python scripts/pssm_pipeline/00_setup.py            # canonicalize + validate query.fasta
python scripts/pssm_pipeline/01_jackhmmer_search.py  # jackhmmer search -> raw MSA
python scripts/pssm_pipeline/02_clean_msa.py          # Stockholm -> N x L matrix, EVEREST A.3.1 filters
python scripts/pssm_pipeline/03_weights.py            # 99%-identity sequence reweighting, Neff
python scripts/pssm_pipeline/04_pssm.py               # per-column aa frequency table
python scripts/pssm_pipeline/05_score.py              # log-odds score every DMS variant
python scripts/pssm_pipeline/06_evaluate.py           # Spearman rho vs DMS, bootstrap CI, scatter plot
```

Pipeline dependency chain (each step consumes only the previous step's
checkpoint, so steps can be re-run individually after changing one of them):
`00 query.fasta` → `01 msa_raw.sto` → `02 msa_clean.npy` → `03 weights.npy` →
`04 pssm.npy` → `05 predictions.csv` → `06 rho + scatter.png`.

Step 1 requires a local UniRef100 FASTA snapshot on disk (see Data acquisition
below) and two constants set by hand at the top of the script: `SEQ_DB`, the
path to the specific year's snapshot to search — there's no year/CLI
parameter, so switching years means editing this constant — and
`BITSCORE_PER_RESIDUE`, configured per protein; thresholds are
length-normalized (bits/residue × query length) rather than e-value-based, so
hit quality stays constant as snapshots grow across years.

## Data acquisition (UniRef100 snapshots)

`scripts/download_uniref100.py` fetches one year's UniRef100 FASTA at a time:

```bash
python scripts/download_uniref100.py --years option2   # this project's locked 13-year set
python scripts/download_uniref100.py --years 2010 2015 --output-dir data/snapshots
```

Historical UniProt releases only ship a combined `uniref{YYYY}_{NN}.tar.gz`
(UniRef50+90+100 together) — no standalone UniRef100 file exists for any
archived year. This script streams that tarball, extracts only the
`uniref100.xml.gz` member (always first in the stream), and aborts the
connection immediately after — UniRef50/90 are never requested.

Two mirrors carry these archives at an identical path layout,
`ftp.uniprot.org` and `ftp.ebi.ac.uk`, and which one is healthy changes without
warning — each has been caught hanging at the TLS handshake, and each has been
seen serving bodies hundreds of times slower than the other. So the script
doesn't hardcode one. At the start of a batch it reads 8 MB from each and pins
whichever is *fastest*, not whichever answers first: a mirror can return
headers promptly while delivering at 0.4 MB/s, which is the difference between
this finishing in an hour and in a fortnight. The chosen mirror is then fixed
for every year in the batch, because a dropped transfer resumes with an HTTP
`Range` request into a still-live decompressor — switching hosts mid-file would
splice two copies together. A mirror dying mid-batch fails those years rather
than rotating; rerun and it re-probes, then resumes from what's on disk.

`--min-free-gb` (default 150) refuses to begin a year's download when the
output volume is below that floor. Peak disk is larger than the final FASTA:
a year's `.xml.gz` is deleted only once its FASTA passes the integrity check,
so concurrent downloads hold several multi-hundred-GB intermediates at once.

Download and parse are two separate, independently-concurrent phases,
connected by the extracted `uniref100.xml.gz` sitting on disk rather than a
live pipe: `--download-workers` (default 4) controls how many years fetch at
once — cheap, network/disk-bound — and `--parse-workers` (default 3) controls
how many `scripts/xml_to_fasta.py` conversions run at once — CPU-bound, and
capped low to keep total memory use predictable. The `.xml.gz` is deleted
once its FASTA is verified. `scripts/xml_to_fasta.py` is a line-oriented
streaming parser (no DOM, memory flat at ~one sequence) and is independently
invocable.

`scripts/run_tier_a_download.sh` and `scripts/run_tier_b_download.sh` launch
the 2010–2018 and 2020/2022/2024/2026 batches respectively as detached
background jobs (`nohup` + PID file under `data/snapshots/`) that survive an
SSH disconnect, logging to `logs/` (gitignored). Check/stop a running download
with the `kill -0`/`kill` commands each script prints on start. Output goes to
`data/snapshots`, which is expected to be a symlink to scratch/NVMe storage —
the root volume does not have room for multi-hundred-GB snapshots.

This repo is checked out on more than one machine (a local Mac and an EC2
instance), and each has its own downloaded snapshots — which years are
present, and where `data/snapshots` actually points, differs per machine and
is deliberately not committed. Check what's on disk (`ls data/snapshots`)
rather than assuming another machine's state.

## Data layout

- `data/protein.fasta` — query: full-length SARS-CoV-2 Spike (1273 aa),
  precursor numbering. Tracked in git (small, canonical input).
- `data/SARS2_RBD_Starr_binding_dms.csv` — Starr 2020 ACE2 binding DMS
  (RBD positions 331–531), same numbering as the query — no coordinate offset
  needed anywhere in the pipeline.
- `data/pssm_pipeline/` — gitignored checkpoints written by pipeline steps
  00–06; regenerate by rerunning the relevant script.
- `data/snapshots/` — gitignored multi-GB UniRef100 FASTA snapshots (one per
  acquired year) plus `.stats.json` sidecars; typically a symlink to
  scratch/NVMe, not committed or kept on the root volume.
- `data/uniprotref_yearly_archive_sizes.csv` — combined UniRef50+90+100
  archive size per year, used for both the growth plot and
  `download_uniref100.py`'s integrity checks (expected cluster counts).
- `calibration.csv` — timed measurements (wall/CPU/RSS/throughput) from
  early conversion and search calibration runs; referenced by the README's
  "Compute breakdown" section, not consumed by any script at run time.

## Key methodology to preserve when modifying the pipeline

- **Bit-score thresholds, not e-values**, for the jackhmmer search: e-value
  cutoffs would silently tighten as later-year snapshots grow, while
  bits/residue × query length keeps hit stringency constant across the
  2010–2026 sweep.
- **Sequence reweighting at 99% identity** (theta=0.01), not the more common
  80%, because two sequences differing by 1% can already have meaningfully
  different fitness for this protein. Neff/L (effective sequences per column)
  is the depth-adequacy check against EVEREST's floor of 1.0.
- **Column/sequence filtering in `02_clean_msa.py`** (drop columns >50% gaps,
  drop sequences <50% query coverage) is computed against the *original*
  query positions independently for both filters, not against each other's
  output — preserve this independence if the filters are ever touched.
- Coordinate mapping between the MSA, PSSM, and DMS is currently offset-free
  by construction (single-query jackhmmer profile ⇒ match column *i* is query
  position *i*; DMS and query already share numbering) — if the query or
  search strategy changes such that this no longer holds, every downstream
  step's coordinate assumptions need re-verification.
