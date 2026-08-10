# `sweep_results.csv` — data dictionary and interpretation notes

One row per UniRef100 snapshot year. Produced by `scripts/sweep/collect.py`,
which re-derives the table from the per-year pipeline checkpoints under
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

## Columns

### Identity / database
| Column | Meaning |
|---|---|
| `tag` | Run directory name under the sweep root. Equals the year. |
| `year` | UniRef100 release year (the January release, `uniref100_<year>_01`). |
| `snapshot_bytes`, `snapshot_gb` | Size of the FASTA actually searched. |
| `db_n_seqs`, `db_n_residues` | Sequence / residue count of that snapshot, from the download step's `.stats.json`. Prefer these over bytes for a "database size" axis. |

### Search (step 01, jackhmmer)
| Column | Meaning |
|---|---|
| `bitscore_per_residue` | Inclusion threshold in bits per residue. **0.3 for every row** — held constant by design. |
| `query_length` | 1273 (full-length Spike, precursor numbering). Constant. |
| `threshold_bits` | `bitscore_per_residue × query_length` = 381.9. Constant. A raw bit score, not an E-value, so the cutoff does not silently tighten as snapshots grow. |
| `jackhmmer_elapsed_s` | Wall-clock seconds. **See the compute-cost warning below — do not treat as a clean benchmark.** |
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
| `clears_depth_floor` | Whether `Neff_over_L ≥ 1.0`, EVEREST's selection threshold. **Expected False throughout — see caveats.** |
| `Neff_at_90pct_identity` | Reliability metric from Methods A.6.1; paper's threshold is 30. |
| `clears_reliability` | Whether that threshold is met. |
| `n_singleton_sequences` | Sequences in a cluster of one — a redundancy indicator. |

### Scoring (step 05)
| Column | Meaning |
|---|---|
| `n_variants` | DMS variants attempted. **3802 every year** — the DMS is fixed. |
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

**The decline is substantially a coverage effect, not a quality effect.** See
caveat 7 — `L_final` is the confound, and `spearman_rho_excl_imputed` behaves
differently enough that the two columns support different conclusions. Do not
report the rho curve alone.

---

## Interpretation caveats

**1. Every snapshot here is pre-pandemic.** SARS-CoV-2 emerged in late 2019; the
available years are 2010–2018. The query sequence therefore appears in *none* of
these databases. What accumulates across the series is other coronaviruses
(SARS-CoV-1, MERS from 2012, bat/civet relatives), not the query's own lineage.
**A flat, weakly-rising, or declining curve is a legitimate scientific result**,
not a broken pipeline — it would say the marginal value of pre-pandemic
sequencing for this specific target is small or negative, which is precisely the
question the project asks. Do not present such a curve as a failure, and do not
stretch axes to manufacture a trend.

The observed curve declines (see Headline result). Read together with caveat 7,
the defensible reading is that growth in these pre-pandemic databases recruits
increasingly distant homologs, which degrades alignment *coverage* of the RBD
faster than it improves the substitution statistics — so the marginal value of
this data for this target is at best zero and plausibly negative.

**2. `spearman_rho` vs `spearman_rho_excl_imputed` measure different things.**
Imputed variants all receive one constant value, which carries zero rank
information, so they dilute the headline rho toward 0. `spearman_rho` is the
honest end-to-end number and is comparable across years (fixed n=3802).
`spearman_rho_excl_imputed` better reflects the model where it actually makes
predictions — but its sample differs year to year, so cross-year comparison of
that column is not strictly like-for-like. In this completed sweep that
caveat is much stronger than "not strictly": `n_excl_imputed` ranges from 2211 to
3175 and the *difficulty* of the retained subset shifts systematically with it
(caveat 7). Plot both; if plotting one, use `spearman_rho` and show the imputed
fraction alongside it.

**3. `Neff_over_L` is expected to sit far below EVEREST's depth floor of 1.0**
(2010 is 0.252). Under EVEREST's own selection heuristic these alignments would
not be chosen for downstream modeling. This is a real limitation on how much
weight the absolute rho values can bear, and should be stated wherever the
numbers are presented.

**4. These do not reproduce EVEREST's published numbers, by construction.** One
dated snapshot instead of the paper's full retrieval, and one bit-score threshold
instead of its eighteen alignment variants. PSSM is also the weakest
alignment-based model in the paper (avg rho ≈ 0.38 vs 0.44 for EVE). A
single-protein rho in the 0.2–0.5 band is plausible; **rho ≈ 0 would signal a
bug, most likely in coordinate mapping.**

**5. `jackhmmer_elapsed_s` is NOT a clean compute benchmark.** These runs were
executed six-at-a-time on a 16-vCPU box. A measured control (the same 2010
database run alone vs. six-way concurrent) showed **223.8 s → ~403 s, a 1.80×
inflation** purely from contention. Worse, concurrency was not even constant
across the sweep: eight years were queued six-at-a-time and the field thinned as
short years finished, so 2018 ran with far less contention than 2015 did. The
column therefore conflates database size with a scheduling history that is not
recorded anywhere. Do not fit a scaling law of runtime vs. database size to it,
and do not quote it as a per-year compute cost. `README.md`'s compute-cost
section has properly isolated measurements for that purpose.

**6. Year coverage stops at 2018 and is not extendable from this machine
as-is.** The 2020 snapshot exists but is 0 bytes (a killed download); 2022/2024/
2026 were never fetched. Any post-pandemic point — which is where the
interesting signal would be, since it is the first time the query's own lineage
enters the database — requires re-running the Tier B download first.

**7. `L_final` confounds the headline curve, and is the most important thing to
handle before plotting.** `L_final` is not stable across the series: it sits at
872–875 for 2012–2015, then falls to 829/818/820 for 2016/2017/2018, and is
844/837 for 2010/2011. Because imputation is driven entirely by which columns
survive, `imputed_frac` tracks it inversely — 16.5–17.0% in the high-`L_final`
years, 33–42% in the low ones. Three consequences:

- Headline rho is depressed in exactly the years with low `L_final`, because
  more variants receive the constant fill that carries no rank information. Much
  of the 2015 → 2018 decline is this, not a change in the model.
- `rho_excl_imputed` moves *opposite* to rho across that boundary: it jumps from
  0.1143 (2015) to 0.1711 (2016). When the gap filter drops ~45 more columns, the
  survivors are the best-covered ones and the retained variant subset gets easier.
  This is a difficulty shift in the denominator, not an improvement.
- 2016 and 2017 have nearly identical `rho_excl_imputed` (0.1711 vs 0.1708) while
  their headline rho differs by 0.012 — the difference is fully explained by
  `imputed_frac` rising 37.4% → 41.9%.

The mechanism is coherent: a larger database recruits more distant homologs
(`N_final` 451 → 2745), which makes the alignment gappier, which drops columns,
which raises imputation. Database growth is degrading *coverage* of the RBD while
per-prediction accuracy stays roughly flat within a coverage regime.

**Practical guidance:** plot `L_final` and `imputed_frac` as companion panels
rather than mentioning them in a caption, and treat any claim about the rho trend
as provisional until it survives conditioning on `L_final`. With only nine points
split unevenly across coverage regimes, database size and column retention are
not cleanly separable in this dataset — say so rather than picking whichever
column tells the tidier story.

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
