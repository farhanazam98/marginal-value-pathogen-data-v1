# `sweep_results.csv` — data dictionary

One row per (protein, UniRef100 snapshot year, DMS assay), keyed by
`(protein, tag, dms_id)`. Steps
00–04 build one PSSM per year; steps 05–06 score every assay named by the
active protein config against it, so a year with two assays produces two rows
sharing all step 00–04 columns. Produced by `scripts/sweep/collect.py`, which
re-derives the table from the per-year pipeline checkpoints under
`/mnt/scratch/sweep/<year>/data/pssm_pipeline/`. Regenerate at any time with:

```bash
conda activate marginal-value-pathogen-data
python scripts/sweep/collect.py
```

The experiment: for each yearly snapshot of UniRef100, run the full PSSM
pipeline (`scripts/pssm_pipeline/00_*.py` … `06_*.py`) against SARS-CoV-2 Spike
and Spearman-correlate the predicted mutation effects against the Starr 2020
ACE2-binding DMS. Every year uses an identical pipeline and an identical
bit-score threshold (0.3 bits/residue); **the snapshot is the only variable.**

---

## Terms

No biology background assumed.

- **Homolog** — a protein in another organism descended from a shared
  ancestor. Here: spike proteins of other coronaviruses.
- **MSA (multiple sequence alignment)** — homologs stacked so equivalent
  positions line up in columns, with gaps inserted where lengths differ.
- **jackhmmer** — the search tool that finds homologs in a snapshot. The
  **bit-score threshold** is the bar for "related enough to include"; a raw
  bit score, held at 0.3/residue for every year so the filter never moves.
- **PSSM (position-specific scoring matrix)** — the model. For each position,
  count how often each amino acid appears there across the aligned homologs.
  Rare substitution at a well-conserved position → predicted damaging.
- **`Neff`** — homolog count corrected for redundancy. Databases hold many
  near-identical copies of heavily-sequenced organisms, so sequences ≥99%
  identical are clustered and each cluster counted once: "how many genuinely
  distinct relatives."
- **`L_final`** — alignment columns surviving the gap filter. Columns >50%
  gaps are too sparse to learn from and get dropped.
- **DMS (deep mutational scanning)** — the ground truth. Starr et al. built
  ~3800 single RBD mutations and measured ACE2 binding for each.
- **Imputation** — filling a missing value with a stand-in. A variant whose
  column was dropped has no prediction, so it receives the mean prediction.
- **Spearman rho** — rank correlation: did the model order the mutations
  correctly, not were its numbers right. 1 = perfect, 0 = no relationship.
- **Bootstrap CI** — resample the 3802 variants 10,000 times, recompute rho
  each time, keep the middle 95%. Wide or overlapping intervals mean two
  years cannot be called different.

---

## Columns

### Identity / database
| Column | Meaning |
|---|---|
| `protein` | Which protein this row's PSSM was built for, from the sandbox path (`<protein>/<year>`). The protein config's filename stem (e.g. `spike`). |
| `tag` | Run directory name under the sweep root. Equals the year. |
| `dms_id` | Which DMS assay this row scores (e.g. `starr_binding`, `starr_expression`), from the protein config. The one column that distinguishes multiple rows of the same year. |
| `year` | UniRef100 release year (the January release, `uniref100_<year>_01`). |
| `snapshot_bytes`, `snapshot_gb` | Size of the FASTA actually searched. |
| `db_n_seqs`, `db_n_residues` | Sequence / residue count of that snapshot, from the download step's `.stats.json`. Prefer these over bytes for a "database size" axis. |

### Search (step 01, jackhmmer)
| Column | Meaning |
|---|---|
| `bitscore_per_residue` | Inclusion threshold in bits per residue. **0.3 for every row** — held constant by design. |
| `query_length` | 1273 (full-length Spike, precursor numbering). Constant. |
| `threshold_bits` | `bitscore_per_residue × query_length` = 381.9. Constant. A raw bit score, not an E-value, so the cutoff does not silently tighten as snapshots grow. |
| `jackhmmer_elapsed_s` | Wall-clock seconds. **Not a clean benchmark** — runs were executed six-at-a-time on a 16-vCPU box (measured contention inflation: 1.80×), and concurrency thinned as short years finished. Do not quote as a per-year compute cost; `README.md` has isolated measurements. |
| `jackhmmer_rounds` | Iterations used (cap is 5). |
| `jackhmmer_converged` | True only if a round added exactly 0 new targets. Commonly false at the 5-round cap while oscillating around ~1 new target; not a failure. |
| `n_hits` | Significant hits found. |
| `n_alignment_rows` | Rows in the raw Stockholm alignment (`n_hits` + the query seed row). |

### Alignment cleaning (step 02)
| Column | Meaning |
|---|---|
| `N_raw`, `L_raw` | Sequences / columns before filtering. `L_raw` is always 1273. |
| `N_final` | Sequences surviving the ≥50% query-coverage filter. |
| `L_final` | Columns surviving the ≤50%-gap filter. **Key driver of everything downstream:** DMS variants at positions whose column was dropped cannot be scored and get imputed. |

### Sequence weighting (step 03)
| Column | Meaning |
|---|---|
| `theta` | 0.01 → cluster at 99% identity. Constant. |
| `Neff` | Effective sequence count: sum of weights, where each sequence's weight is 1/(size of its ≥99%-identity cluster). Corrects for databases oversampling intensively-sequenced lineages. |
| `Neff_over_L` | `Neff / L_final`. The standard alignment-depth statistic. |
| `clears_depth_floor` | Whether `Neff_over_L ≥ 1.0`, EVEREST's selection threshold. **False for every year** — under EVEREST's own heuristic these alignments would not be selected for downstream modeling, which limits how much weight the absolute rho values can bear. |
| `Neff_at_90pct_identity` | Reliability metric from Methods A.6.1; paper's threshold is 30. |
| `clears_reliability` | Whether that threshold is met. |
| `n_singleton_sequences` | Sequences in a cluster of one — a redundancy indicator. |

### Scoring (step 05)
| Column | Meaning |
|---|---|
| `n_variants` | DMS variants attempted. Fixed per assay across years (Starr binding 3802, expression 3798); it varies only between assays, not between snapshots. |
| `n_scored_directly` | Variants landing on a surviving MSA column, i.e. genuinely predicted. |
| `n_imputed`, `imputed_frac` | Variants whose column was dropped in step 02. These receive a **constant** fill value, not a prediction. |
| `imputed_value` | The constant used (the mean of the directly-scored predictions). |
| `wt_wt_all_zero` | Sanity flag: wild-type→wild-type must score exactly 0. Must be True. |
| `predicted_score_mean`, `predicted_score_std` | Distribution of predicted scores. |

### Evaluation (step 06)
| Column | Meaning |
|---|---|
| `n_joined` | Variants matched between predictions and DMS. Should equal 3802. |
| `n_dropped_from_dms` | Join losses. Should be 0. |
| `spearman_rho` | **Headline metric.** Spearman correlation over *all* 3802 variants, imputed included. |
| `spearman_pvalue` | p-value for that correlation. |
| `bootstrap_ci_95_lo`, `bootstrap_ci_95_hi` | 95% CI from 10,000 resamples (seed 0). Use for error bars. |
| `spearman_rho_excl_imputed` | Spearman over only the directly-scored variants. |
| `n_excl_imputed` | Sample size for that number (= `n_scored_directly`). |
| `status` | `DONE`, or `FAILED:<step>` for an incomplete run. **Filter to `DONE` before plotting.** |

---

## Headline result (all 9 years complete, 2026-08-10)

| year | GB | N_final | L_final | Neff/L | imp% | rho | 95% CI | rho excl. imp |
|---|---|---|---|---|---|---|---|---|
| 2010 | 4.1 | 451 | 844 | 0.252 | 33.0 | 0.1750 | [0.1417, 0.2083] | 0.2031 |
| 2011 | 4.9 | 494 | 837 | 0.268 | 39.4 | 0.1590 | [0.1265, 0.1917] | 0.2013 |
| 2012 | 6.6 | 606 | 875 | 0.339 | 17.0 | 0.1542 | [0.1226, 0.1863] | 0.1681 |
| 2013 | 8.6 | 767 | 872 | 0.376 | 16.5 | 0.1480 | [0.1163, 0.1796] | 0.1568 |
| 2014 | 14.4 | 1000 | 872 | 0.449 | 16.5 | 0.1364 | [0.1052, 0.1691] | 0.1493 |
| 2015 | 21.9 | 1287 | 873 | 0.518 | 16.5 | 0.1047 | [0.0733, 0.1374] | 0.1143 |
| 2016 | 31.9 | 1665 | 829 | 0.604 | 37.4 | 0.1285 | [0.0954, 0.1617] | 0.1711 |
| 2017 | 41.9 | 2164 | 818 | 0.717 | 41.9 | 0.1162 | [0.0820, 0.1495] | 0.1708 |
| 2018 | 58.8 | 2745 | 820 | 0.893 | 41.3 | 0.0996 | [0.0649, 0.1335] | 0.1410 |

**Headline rho declines as the database grows** — 0.175 at 4.1 GB to 0.100 at
58.8 GB, a 43% relative drop, while alignment depth more than triples
(`Neff_over_L` 0.252 → 0.893). The 2010 and 2018 CIs are disjoint. The decline is
not monotone: 2016 breaks it upward.

`L_final` is not stable across the series (872–875 for 2012–2015; 818–844 for the
other years) and `imputed_frac` tracks it inversely, so database size and column
retention are not cleanly separable in these nine points.

---

## Provenance

- The 2010 row was promoted from a six-way concurrency probe in which all six
  independent runs produced **bit-identical** results (`rho =
  0.17495952679056165`), matching the repository's pre-existing 2010 checkpoint
  and the separately-run Mac result. This validates both the pipeline's
  determinism and the sandbox isolation used for the sweep.
- No tracked pipeline script was modified to produce this sweep. Each year ran in
  an isolated sandbox whose `data/uniref100_2010.fasta` symlink pointed at that
  year's snapshot — the filename is retained because `SEQ_DB` in
  `01_jackhmmer_search.py` hardcodes it.
- All nine runs completed 2026-08-10 (driver finished 11:28 UTC). Every row was
  verified to have `n_variants` = `n_joined` = 3802, `n_dropped_from_dms` = 0,
  `wt_wt_all_zero` = True, `threshold_bits` = 381.9, and `jackhmmer_rounds` = 5.
  `jackhmmer_converged` is False for 2010, 2011, and 2013 — the 5-round cap
  reached while oscillating around ~1 new target, which is expected and not a
  failure (see the column notes).
