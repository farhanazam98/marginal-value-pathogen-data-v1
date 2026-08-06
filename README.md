# marginal-value-pathogen-data-v1

**Objective:** Measuring whether removing sequences >90% identical to the SARS-CoV-2 RBD changes PSSM's ability to predict the Starr 2020 ACE2 binding assay, across cumulative database snapshots from 2010 to 2025, with a depth-matched control.

## Setup

Analysis scripts depend on `matplotlib`, distributed via conda. On a fresh machine
(e.g. a new EC2 instance):

```bash
./setup.sh
conda activate marginal-value-pathogen-data
```

`setup.sh` installs Miniconda if it isn't already present, accepts the conda channel
Terms of Service (required non-interactively), and creates/updates the
`marginal-value-pathogen-data` environment from `environment.yml`. Re-running it is safe
and just updates the existing environment.

## Milestones:
### Monday
Come up with objective. Find and download, and the Starr 2020 DMS data for SARS-COV-2 in the [priority-viruses](https://github.com/debbiemarkslab/priority-viruses) repo and parse the N501Y mutation. 

*Deliverable: CSV file with DMS data from Starr 2020 experiment for SARS-COV-2 under data/ and a script that parses the N501Y mutation*

### Tuesday
Find databases needed for SARS-COV-2 multiple sequence alignment (MSA), segmented by year. 

*Deliverable: Database snapshots located, along with the size for each respective database. This will give us an idea of how long it will take to download all the databases on an EC2 instance.* 

As a proxy for overall database growth (and download time), here's the size of
UniProt's combined UniRef50+90+100 [archive](https://ftp.uniprot.org/pub/databases/uniprot/previous_major_releases/release-2011_01/uniref/) across the first release of each year from 
2010–2026 (`see data/uniprotref_yearly_archive_sizes.csv`). Note that this is the size for the compressed file, which may be ~7-8x as large when uncompressed. 

![UniProt UniRef archive growth, 2010-2026](data/uniprotref_archive_growth.png)

### Wednesday
For SARS-COV-2, build MSA using JackHMMER, fit a PSSM, and compare against DMS, calculating Spearman Coefficient

*Deliverable: A score for the mutation effect prediction of a single protein, and spearman coefficient. *

## Pipeline (alignment-based zero-shot variant effect prediction)

Minimal reimplementation of the alignment-based half of the EVEREST pipeline
(Gurev/Youssef/Marks, bioRxiv 2025.08.04.668549) for SARS-CoV-2 Spike: retrieve
homologs → build MSA → reweight sequences → fit a PSSM → score all single
substitutions in the RBD → Spearman-correlate against the Starr 2020 DMS.

Pipeline scripts live under `scripts/pssm_pipeline/` (numbered `00_*.py` ...
`06_*.py`, one per step, runnable standalone from the repo root). Intermediate
checkpoints live under `data/pssm_pipeline/`, so each step's output can be
inspected independently without rerunning anything.

Query: `data/protein.fasta`, full-length Spike (1273 aa), precursor numbering
(starts at the initiator Met). DMS: `data/SARS2_RBD_Starr_binding_dms.csv`
(Starr 2020 ACE2 binding assay, positions 331-531), which uses this same
numbering — no coordinate offset needed between query and DMS.

Caveats to keep in mind throughout: the search runs against a single dated
UniRef100 snapshot rather than the paper's full pipeline, and only one
bit-score threshold is used here instead of the paper's eighteen alignment
variants, so this will not reproduce EVEREST's reported numbers. PSSM is also
the weakest of the alignment-based models in the paper (avg rho ~0.38 vs 0.44
for EVE) — a single-protein rho in the 0.2-0.5 range is plausible; rho near 0
signals a bug, most likely in coordinate mapping.

### Step 0 — `scripts/pssm_pipeline/00_setup.py`

Canonicalizes the query sequence into `data/pssm_pipeline/query.fasta` and
validates it before any downstream step depends on it: right length, standard
amino acids only, correct start/end.

Run: `python scripts/pssm_pipeline/00_setup.py`

Results:
- L = 1273 (matches full-length Spike)
- Starts `MFVFLVLL...`, ends `...VKLHYT` — correct precursor boundaries
- No non-standard amino acid characters
- Residue composition looks like a typical viral glycoprotein (elevated S/T/L/N)

All sanity checks passed.

### Step 1 — `scripts/pssm_pipeline/01_jackhmmer_search.py`

Runs `jackhmmer` (HMMER suite) locally against a local UniRef100 FASTA
snapshot (`data/uniref100_2010.fasta`), which iteratively searches the
database for homologs (sequences descended from a common ancestor). Each
iteration builds a profile from the hits found so far and re-searches with it,
so distantly related sequences that a single-pass search would miss can still
be detected. The inclusion threshold controls how distant a homolog can be
before it's excluded; EVEREST length-normalizes it (bits/residue x L) so long
and short proteins get comparably strict cutoffs. This step fixes the homolog
set that every later step depends on.

This replaces the earlier EBI HMMER web-service version. The web service runs
the same `jackhmmer` binary, but its databases have no UniRef100 option, so it
fell back to `uniprot`; running locally lets us search the actual UniRef100
release the paper uses, and collapses the whole submit/poll/download HTTP dance
into a single subprocess call (jackhmmer writes the alignment straight to disk
with `-A`). The four API bit-score parameters map directly to CLI flags:
`T`/`incT`/`domT`/`incdomT` -> `-T`/`--incT`/`--domT`/`--incdomT`. Requires the
HMMER suite on PATH (`conda activate marginal-value-pathogen-data`).

Run: `python scripts/pssm_pipeline/01_jackhmmer_search.py`

Results (0.3 bits/residue, UniRef100 2010 snapshot):
- Threshold: 0.3 bits/residue x L=1273 = 381.9 bits (bitscore mode, not
  e-value mode: a bit score is database-size-independent, whereas an e-value
  cutoff would silently tighten as the 2010->2026 snapshots grow -- so a fixed
  bits/residue keeps hit quality constant across snapshots). Configure
  `BITSCORE_PER_RESIDUE` per protein.
- Effectively stable but not formally converged within the 5-iteration cap:
  new targets per round were 322 -> 216 -> 1 -> 12 -> 1. jackhmmer only prints
  CONVERGED when a round adds exactly 0 new targets, so it hits the -N 5 cap
  while oscillating around 1; a slightly higher -N would likely tip it to a
  formal convergence.
- 549 significant hits; 553 rows in the Stockholm alignment.
- Top 10 hits: bit scores ~1835-1837, e-value 0, all annotated "Spike
  glycoprotein" (SARS-CoV-1 and related coronaviruses -- the pre-pandemic
  homologs available in 2010).
- No verbatim database match for the query -- expected, since SARS-CoV-2
  post-dates the 2010 snapshot. The one query-identical alignment row is the
  jackhmmer *seed* itself (named `my_protein`, not a `UniRef100_` accession):
  jackhmmer's output MSA is by construction "the query aligned with its
  homologs", so the query is always present as a row regardless of whether a
  matching database entry exists. Step 2's find-the-query-row step relies only
  on this seed row.

Wrote `data/pssm_pipeline/msa_raw.sto` (alignment, consumed by Step 2),
`data/pssm_pipeline/msa_raw_hits.tbl` (per-sequence hit table), and
`data/pssm_pipeline/msa_raw_run_meta.json` (run summary / debug info).

### Step 2 — `scripts/pssm_pipeline/02_clean_msa.py`

Turns the raw Stockholm alignment into a plain N x L positional matrix that
later steps can index into directly. Coordinate mapping is simple here: this
jackhmmer profile was built from a single query sequence (iteration 1), so it
has exactly L=1273 match columns, one per query residue, and match-column
index i *is* query position i with no offset arithmetic needed -- verified by
extracting only the match-column characters from the self-hit row and getting
back the query sequence exactly. Insert columns (residues some hits have that
the query doesn't) are simply discarded.

Then EVEREST's Methods A.3.1 filters apply: columns where most homologs are
gaps get dropped (>50%), since a PSSM position built mostly from placeholders
isn't measuring real conservation; sequences that barely overlap the query
get dropped (<50% coverage), since they contribute little signal but can
still skew frequency counts at the columns they do touch. Both filters are
computed against the original 1273 query positions independently of each
other (not against each other's output), so "50% gaps" and "50% coverage"
mean exactly what they say regardless of filter order.

Run: `python scripts/pssm_pipeline/02_clean_msa.py`

Results:
- Columns: 1273 -> 844 (429 dropped for >50% gaps)
- Sequences: 553 -> 451 (102 dropped for <50% query coverage)
- Final matrix shape: (451, 844)
- Query row: 0 gaps, matches `data/query.fasta` exactly at every one of the
  844 kept positions
- 429 of 1273 columns dropped means many RBD positions (331-531) don't
  survive, so Step 5 imputes a large share of DMS variants (1254/3802)
- Column gap-fraction distribution (p50=0.16, p75=0.55) is bimodal: most
  columns are well-covered, with a second cluster above ~0.55 at the
  worst-covered columns (alignment ends / variable loops)

All sanity checks passed. Wrote `data/pssm_pipeline/msa_clean.npy` and
`data/pssm_pipeline/msa_clean_meta.json`.

### Step 3 — `scripts/pssm_pipeline/03_weights.py`

Raw sequence counts overstate the evidence a database provides, because
databases sample what people happened to sequence, not what evolution
actually explored -- a handful of intensively-surveilled lineages can flood
an alignment with thousands of near-duplicate sequences, drowning out
genuinely distinct ones. Sequence weighting corrects for this: each
sequence's weight is 1 over the size of its near-identical cluster, so ten
copies of the same strain collectively count as one independent observation.
EVEREST clusters at 99% identity (theta=0.01) rather than the more common
80%, because two sequences differing by just 1% can already have meaningfully
different fitness for a protein like this. Neff (sum of weights) is the
"effective" number of independent sequences, and Neff/L is the standard
rule-of-thumb for whether there's enough independent evidence per model
position to fit a reliable per-column frequency estimate.

Run: `python scripts/pssm_pipeline/03_weights.py`

Results:
- N=451, Neff=212.66 at theta=0.01 (99% identity clustering) -- ~47% of the
  raw sequences count as independent evidence once near-duplicates collapse
- **Neff/L = 0.252, below EVEREST's depth floor of 1.0.** Under EVEREST's own
  selection heuristic this alignment would not be selected for downstream
  modeling. Not a bug: unlike the old `uniprot` run (flooded with SARS-CoV-2
  near-duplicates), this alignment is simply *shallow* -- the 2010 snapshot
  predates the virus, so only 451 distant coronavirus homologs are found. The
  5 largest clusters are just 53/52/51/49/46 and 157/451 sequences are
  singletons, so redundancy is modest; there just aren't many sequences.
- Neff @ 90% identity (Methods A.6.1 reliability metric) = 58.56, which
  **clears** the paper's threshold of 30
- The depth-floor failure is a real caveat for interpreting the final
  Spearman correlation, not just a formality

All sanity checks passed. Wrote `data/pssm_pipeline/weights.npy` and
`data/pssm_pipeline/weights_meta.json`.

### Step 4 — `scripts/pssm_pipeline/04_pssm.py`

A PSSM (position-specific scoring matrix) is a per-column amino-acid
frequency table: for each of the 844 surviving alignment columns, how often
does each of the 20 amino acids show up, once each sequence's vote is scaled
by its Step-3 weight so that a cluster of near-duplicates doesn't get to vote
once per member. This is the "site-independent" approximation: it captures
how conserved each position is on its own, but says nothing about
correlations between positions (e.g. compensating mutations), which is
exactly what an actual maximum-entropy model like EVE/Hopf et al. 2017 adds
on top.

Run: `python scripts/pssm_pipeline/04_pssm.py`

Results:
- Every column's 20 frequencies sum to 1 (verified)
- Query row still gap-free (0 gaps), as it must be after Step 2
- WT (wild-type) residue -- the amino acid present at that position in our
  query, as opposed to a mutant substitution -- is the single most frequent
  amino acid in only 310/844 columns (36.7%), far below the old `uniprot`
  run's 95%. The 2010 homologs are distant coronaviruses whose consensus
  often differs from SARS-CoV-2's residue; WT reconciliation still passes
  100% (Step 5), so this reflects real cross-species divergence, not a bug.
- Per-column entropy ranges from 0.76 to 3.72 bits (max = log2(20) = 4.32),
  mean 2.31 bits -- higher than the old run, i.e. a more diverse homolog set
- The 10 most-conserved columns are cysteines (e.g. positions 538, 590, 738,
  749) at freq ~0.92 -- structurally load-bearing disulfide-bonded residues
- The 10 least-conserved columns cluster in the RBD (e.g. positions 436, 460,
  489, 492, 501, 502, 506, 507) -- the most variable region across
  coronaviruses, not noise

All sanity checks passed. Wrote `data/pssm_pipeline/pssm.npy` and
`data/pssm_pipeline/pssm_meta.json`.

### Step 5 — `scripts/pssm_pipeline/05_score.py`

For every variant in the DMS (deep mutational scan -- an experiment that
measures a fitness-like readout for thousands of point mutations in
parallel), the predicted score is `log f(mut, pos) - log f(wt, pos)` using
the Step 4 PSSM: the log-ratio of how often the mutant amino acid appears at
that alignment column versus how often the wild-type does. This treats the
alignment as a record of an experiment natural selection already ran --
amino acids tolerated at a position accumulate there in rough proportion to
how little they hurt fitness, so a higher relative frequency is a cheap
proxy for a milder effect. 

Run: `python scripts/pssm_pipeline/05_score.py`

Results:
- WT-residue reconciliation: 3802/3802 (100%) match -- the DMS's `mutant`
  column and `query.fasta` already share full-spike, 1-indexed coordinates,
  so no offset correction was needed
- 2548/3802 variants scored directly; 1254/3802 (33%) fell outside the 844
  surviving MSA columns and were imputed with the mean predicted score across
  scored variants (-1.69). The large imputed share is a direct consequence of
  the shallow 2010 alignment dropping many RBD columns in Step 2.
- WT->WT sanity check passed: all 135 covered positions score exactly 0.0
- Score distribution: mean -1.69, median -1.92, min -4.85, max 3.73; 80.7% of
  scored variants negative (deleterious-leaning)
- The 10 most-deleterious predictions all sit at position 523 (a conserved
  threonine), tied at -4.85 where the mutant amino acid's weighted count is 0
- The 10 least-deleterious predictions (e.g. A475G at +3.73) are variable,
  mutation-tolerant positions, as expected

All sanity checks passed. Wrote `data/pssm_pipeline/predictions.csv` and
`data/pssm_pipeline/predictions_meta.json`.

### Step 6 — `scripts/pssm_pipeline/06_evaluate.py`

Spearman's rho is Pearson correlation computed on *ranks* rather than raw
values: it asks whether the ordering by predicted score and the ordering by
DMS score move together, ignoring the actual scale of either. That's the
right tool here since the log-odds predicted scores and the DMS's binding
measurement aren't on comparable units -- only their relative order is meant
to carry information. Ties (e.g. Step 5's pseudocount-floor ties) get the
average of the ranks they span, scipy's default and the standard convention.

Run: `python scripts/pssm_pipeline/06_evaluate.py` (needs the project conda
env -- scipy isn't on the system `python3`)

Results:
- Join on (position, wt_aa, mut_aa): 3802/3802 variants matched, 0 dropped
  in either direction -- expected, since `predictions.csv` was built from
  this same DMS file in Step 5, so this join is mostly a consistency check
- **Spearman rho = 0.175** (p=1.4e-27), 95% bootstrap CI [0.142, 0.209] over
  10,000 resamples -- real signal (far from 0), at the low edge of the
  0.2-0.5 plausible range
- rho excluding the 1254 imputed variants: 0.204 (n=2548) -- higher than the
  all-variant number, i.e. the large imputed block dilutes the headline rho
  (unlike the old run's 95 imputed, which barely moved it)
- The scatter plot (`data/pssm_pipeline/scatter.png`) shows the expected
  wedge shape: variants with a high (near-zero) predicted score cluster near
  DMS score 0 (tolerated), while variants at the pseudocount floor span nearly
  the *entire* DMS range -- absence from a shallow 2010 alignment is a weak
  signal on its own, consistent with PSSM being the least informative of the
  alignment-based models

All sanity checks passed. Wrote `data/pssm_pipeline/scatter.png` and
`data/pssm_pipeline/evaluate_meta.json`.

Pipeline steps 0-6 now run end-to-end for one bit-score threshold (0.3
bits/residue). Per the deferred TODO, the next step is re-running Step 1
across the full threshold sweep and applying EVEREST's alignment-selection
heuristic instead of this single hardcoded threshold.

### Thursday
For SARS-COV-2, define the sweeps that we want to run to illustrate the relationship between performance on mutation effect prediction and time (via database snapshots). This should include the compute cost, both in terms of space. 

*Deliverable: A CSV file with runs defined, and the compute cost in terms of CPU hours + RAM for each run*

## Compute breakdown: calibration results

Two measured calibration points now anchor the whole cost model: XML-to-FASTA conversion (Phase 1) and a single `jackhmmer` iteration (Phase 2), both run against the 2010 UniRef100 release (9,808,438 sequences, 3.414 Gaa, 4.08 GB FASTA).

### Headline

**Search is not the bottleneck. Data acquisition and conversion are.** The full six-protein, seventeen-year grid is roughly 28 job-hours of search, about 56 with a 2x margin for later-iteration cost growth. Against that, downloading the archive is 840 GB at minimum and converting it is 15.9 single-threaded core-hours. The compute ask for this project is small.

### Measured

| Quantity | Value |
|---|---|
| Conversion, 2010 | 180.8 s, 99.4 percent single-core, 12 MB peak RSS |
| Sequence count check | 9,808,438 vs 9,808,438 published, exact |
| Search, one iteration, 2010 | 10.59 s wall at `--cpu 4` or `--cpu 10` |
| Search throughput | 322 Maa/s, equivalently 385 MB/s of FASTA per job |
| Peak RSS, search | 30 to 145 MB, scaling with `--cpu` buffers, not database size |
| Threshold sensitivity | 0.1 bits/res 10.59 s (246 hits) vs 0.5 bits/res 10.85 s (81 hits) |

### The four findings that matter

**`jackhmmer` does not scale past about two cores on this workload.** Requesting 1, 4, or 10 threads all consumed 21.6 to 21.9 CPU-seconds and kept only 1.8 to 2.3 cores busy. The per-job rate of 385 MB/s is close to a single-threaded FASTA parse ceiling, which suggests the master thread reading and parsing the database starves the workers. Consequence for instance selection: parallelism must come from running many independent (protein, year) jobs concurrently, not from giving one job a large `--cpu`. There are 102 such jobs and they are fully independent.

**Per-protein threshold changes are free.** Raising the threshold fivefold changed wall clock by 2.5 percent, inside run-to-run noise, while cutting retained hits from 246 to 81. The full database scan happens regardless of threshold, so tuning the threshold per protein carries no compute cost.

**Memory is a non-issue at every stage.** Conversion peaked at 12 MB against a 14 GB input, search at 145 MB. Neither scales with database size. RAM should be selected for page-caching the database, not for the processes.

**Storage throughput, not cores, is the EC2 risk.** The 385 MB/s per-job rate exceeds default gp3 EBS baseline of 125 MB/s by roughly 3x, and six concurrent jobs would demand about 2.3 GB/s aggregate. The Mac benchmark hid this entirely because the 4.08 GB database sits in page cache. On EC2 the year's FASTA must be page-cached or on local NVMe, or search becomes I/O-bound and roughly 3x slower than measured.

### Projections to the full study

Derived constants: 0.4159 GB FASTA per million sequences, 18.43 s conversion per million sequences, compressed-XML to FASTA byte ratio 1.53, search throughput 322 Maa/s per job.

| Quantity | Projection |
|---|---|
| Total sequences, 17 archived years | 3,096 M |
| Total residues | 1,078 Gaa |
| Total FASTA if all years retained | 1.29 TB |
| Compressed download, uniref100 only | 840 GB |
| Conversion, all years, single-threaded | 15.9 h (parallel across years: about 2.4 h) |
| Search, 6 proteins x 17 years x 5 iterations | 28 job-hours, 56 with margin |
| Largest single year (2026) | 180 to 198 GB FASTA, 0.78 h search per protein |

Size model validated independently against UR100P at 0.4082 GB per million sequences, a 1.9 percent deviation.

### Open risks

The archive publishes combined UniRef tarballs rather than standalone `uniref100.xml.gz` for historical years, on both mirrors (confirmed in Phase 3). UniRef100 is a measured 43.5 percent of the tarball, but it is always the first member in the tar stream, so a streaming extract that aborts once it clears that member avoids downloading the rest. See [Phase 3: acquisition study](#phase-3-acquisition-study) below — the 840 GB figure stands, the 1.93 TB worst case does not happen.

The two-core ceiling was measured on Apple Silicon with a bioconda arm64 build. Whether it reproduces on EC2 should be confirmed with one short validation run before the instance type is finalised.

## Phase 3: acquisition study

Phase 1/2 established that search is cheap and acquisition is the real cost driver. Phase 3 resolved the two open acquisition risks flagged above (standalone-file availability, tarball size) and locked in the year set and retention policy for Phase 4.

### Headline

**The 1.93 TB worst case does not happen.** UniRef100 is never published standalone in an archived release, on either mirror, but it is always the first member inside the combined tarball, so a streaming extract that aborts once it clears that member transfers only the uniref100 share of the archive. Total download across all 17 archived years: **829 GB**, matching the original 840 GB estimate. **Year set: Option 2** (13 years, dense 2010–2018, sparse 2020–2026, 21.3 h download). **Retention: S3 Standard**, $29/month for the full corpus, chosen over local EBS ($102/month for identical data) and over deleting outright (free, but hours to re-pull a year from the public mirror versus an estimated few minutes from same-region S3).

### Measured (2026-07-31, live against both mirrors)

| Quantity | Value |
|---|---|
| Standalone `uniref100.xml.gz`, any archived release | Not available, either mirror (checked 2011, 2018, 2025/2026) |
| uniref100 member position in combined tar | First (confirmed: 2011 tarball, tar header read from a truncated stream) |
| uniref100 fraction of combined tarball | 43.47 percent (2011, byte-exact), consistent with prior 43.5 percent (2010/2011) |
| Download throughput, UniProt, sustained | 6.45 MB/s, single connection |
| Download throughput, EBI, sustained | 6.18 MB/s, single connection |
| Aggregate throughput, 3 parallel connections | 7.41 MB/s (+15 percent, not +3x — bandwidth-capped, not server-throttled) |
| Combined-tarball total, all 17 years, byte-exact via HTTP HEAD | 1,829.8 GB (supersedes the earlier back-calculated 1.93 TB) |

### Decisions

**Year set: Option 2 (13 years)** — 2010–2018 every year, then 2020, 2022, 2024, 2026. 21.3 h download / 495 GB, versus 35.7 h / 829 GB for the full 17-year grid. Every included-year gap is capped at 2 years. A cheaper 10-year skeleton considered alongside it had a 3-year gap sitting exactly where the corpus is largest and a saturation point is most plausible, which was reason enough to reject it.

**Retention: S3 Standard for all acquired years' FASTA.** $29/month for the full 1.27 TB corpus. Rejected local EBS retention ($102/month for identical data, and it ties the corpus to a kept-around instance/volume) and rejected delete-after-search (free, but a protein added after the fact means re-pulling affected years from the public mirror at 6.45 MB/s — hours per year, versus an estimated few minutes per year to re-pull the same data from same-region S3). Given the marginal monthly cost, retention is worth buying as insurance against the protein set changing.

**Protein set locked by 2026-08-03**, ahead of Phase 4 execution. This is a hard constraint of the acquire-once/search-everything/then-decide economics: once a year's FASTA is deleted, adding a protein later means paying the re-acquisition cost again for every year that wasn't retained.

### Open risk carried into Phase 4

**Download throughput is unverified for EC2.** The 6.45 MB/s figure, and the fact that parallel connections didn't scale (7.41 MB/s aggregate at 3x), were both measured from a sandboxed dev environment, not EC2, and look like a client-side bandwidth ceiling rather than a mirror-side limit. Every download-hours figure in Phase 3 is downstream of this number and should be re-measured on the actual instance before the Phase 4 schedule is finalised — it could move by an order of magnitude in either direction. This is the one number in Phase 3 worth re-checking before trusting the schedule; nothing else here is expected to move it materially.

(S3-to-EC2 retrieval speed and the exact AWS region are still open, but low-stakes: both get settled as a byproduct of picking an instance type in Phase 4, and neither is expected to change the retention or year-set decisions above.)

Full per-year, per-option, and per-retention-tier numbers, each tagged with an evidence source (`measured`/`listed`/`extrapolated`/`assumed`), are in `storage_costs.csv`.

## Phase 4: instance, cost, and schedule

Feasibility is settled. This phase picks hardware, prices the project, and puts a schedule with off-ramps against it. The recommendation, cost, and top risks are summarized below.

### Recommendation

Run one `m8gd.4xlarge` in us-east-1, 16 vCPU, 64 GiB, 950 GB NVMe instance store, $0.923 per hour on demand. Acquisition and conversion fuse into a single stream: pull the tarball, abort after the uniref100 member, decompress and parse in flight, write only FASTA. Compressed XML never touches disk. That drops peak local storage from 408 GB to 250 GB, and it makes conversion disappear from the schedule entirely, because one parse worker sustains 14.73 MB/s of compressed input against a 6.45 MB/s arrival rate and so sits 44 percent busy on a single core. The rule if bandwidth turns out higher is `workers = ceil(observed_MBps / 14.73)`, with the caveat that six concurrent search jobs already consume 12.7 of the 16 vCPUs, leaving parse headroom for only about 29 MB/s; above that, acquisition and search should be staggered rather than overlapped.

Local NVMe rather than more RAM is the position. A single search job streams FASTA at 385 MB/s and six want 2.3 GB/s, which exceeds gp3's 1 GB/s ceiling; provisioning gp3 to that ceiling costs $52 per month, more than the entire instance bill, while `m8gd.4xlarge` bundles the NVMe free. The memory-heavy alternative (`r8g.8xlarge`, 256 GiB, enough to page-cache the largest year) costs $1.89 per hour to solve a problem that is worth hours, not days, so it is not worth buying. Zero GPUs: HMMER's kernel is hand-written CPU SIMD with no GPU path, and nothing else in the pipeline is a neural network.

### Cost

$58 in month one, then $17 per month to retain. That is 39.3 instance-hours at $0.923 ($36), a 50 GB gp3 root ($4), and S3 Standard on 758 GB of FASTA ($17 per month). With a 2x allowance for re-runs and first-time debugging, budget $120. Minimum viable is $26: `c8g.2xlarge`, the 10-year skeleton, gp3, delete after searching. The frontier is flat in protein count, because acquisition does not depend on the protein set and search wall time stays at 2.7 hours whether one protein runs alone or six run concurrently. Within any plausible budget the protein count never constrains the year count, and even the full 17-year grid comes in at $83.

### Schedule

Five weeks, symposium Tuesday 2026-09-01, poster time protected. The schedule is built as a thin vertical slice first rather than stage by stage, so that a complete scientific result exists early and everything after it is breadth rather than validity.

Acquisition splits into two tiers because cost is so unevenly distributed. Tier A is 2010 through 2018, which is 120 GB and 5.2 hours of download, 24 percent of the budget, and yields 9 of the 13 curve points. Tier B is 2020 through 2026, which is 375 GB and 16.1 hours for the remaining 4 points. Week 5 provisions the instance, runs the validation protocol, acquires Tier A, and produces a 9-point rho-versus-year curve for Spike across all three arms. Week 6 acquires Tier B and completes the 13-point single-protein curve, which is the first presentable result and also the fallback for everything downstream. Week 7 scales to 6 proteins and computes the full rho surface. Week 8 is cross-protein analysis and a poster draft frozen by 2026-08-28. Week 9 is print and present. A 2.5x inexperience multiplier is applied to every human task and 1.0x to machine time, since bandwidth does not care who is driving.

The protein set no longer has to hard-lock before acquisition. That constraint applied to the delete-after-search policy; S3 retention was bought instead, so adding a protein later re-reads FASTA from S3 in minutes rather than re-pulling from the public mirror in hours. Only the lead protein is needed on day one, and the full set can be chosen in week 6 once the analysis is understood on one protein.

The schedule does not branch meaningfully on bandwidth, which is the surprising result. At the sandbox rate acquisition is 21.3 hours, one unattended overnight run. At a plausible 50 MB/s it is 2.8 hours, and at 200 MB/s it is 42 minutes, at which point parse and then search become the binding stages instead. No branch requires cutting scope. Scope cuts are triggered by debugging overrun, not by the network, and the fallback ladder runs cheapest first: drop 2026, then 2024, then fall back to the 10-year skeleton, then reduce the protein count, then a cross-shaped design of one protein across all years plus one year across all proteins.

### Top three risks

Download bandwidth on EC2 is still unmeasured and remains the dominant uncertainty, spanning a 31x plausible range. It is also nearly harmless: across that entire range plus a 0.33x to 3x band on search throughput, total pipeline wall clock moves between 21.3 hours and 0.9 hours, and neither endpoint threatens a five-week schedule. Step 1 of the validation protocol collapses it in five minutes, and search throughput only starts to matter above 50 MB/s, where it overtakes download as the binding stage.

The two-core search ceiling was measured on Apple Silicon and the concurrency plan assumes it transfers to Graviton4. Step 2 of the protocol tests it directly; if `--cpu 8` delivers more than about 2.5 effective cores, the right shape is fewer and fatter jobs and the vCPU budget needs recomputing before bulk acquisition.

Inexperience is the real schedule risk. Machine time is roughly a day and the human build is roughly four weeks, so every slip is a human slip. The gate that matters is week 5: if the 2010 reproduction has not passed by Friday 2026-08-07, start the fallback ladder rather than compressing weeks 7 and 8.

