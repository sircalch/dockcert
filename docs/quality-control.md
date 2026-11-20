# How the checks work

Each check returns a status (`PASS`, `WARNING` or `FAIL`) with a short diagnostic. The overall status is the worst of the individual statuses. A `WARNING` asks you to look at a result; a `FAIL` means the evidence does not support the claim as it stands.

## Enrichment metrics

For a table of scores and active/decoy labels, DockCert computes the metrics below, each judged against a pass and a warning threshold (defaults; they are arguments of [`assess_docking_quality`](api.md)).

| Metric | What it measures | `PASS` if | `WARNING` if |
|---|---|---|---|
| ROC-AUC | Probability that an active is ranked above a decoy | \(\ge 0.70\) | \(\ge 0.60\) |
| BEDROC (\(\alpha = 20\)) | Early recognition, weighting the top of the ranked list | \(\ge 0.50\) | \(\ge 0.30\) |
| EF 1 % | Hit rate in the top 1 % relative to random selection | \(\ge 5\) | \(\ge 2\) |
| EF 5 % | Same, for the top 5 % | \(\ge 3\) | \(\ge 1.5\) |
| Optimal MCC | Best Matthews correlation coefficient over score cut-offs | \(\ge 0.40\) | \(\ge 0.20\) |

Below the warning threshold the status is `FAIL`. The precision–recall AUC, RIE and log AUC are also computed and reported. Functions: [`calculate_roc_auc`](api.md), [`calculate_pr_auc`](api.md), [`calculate_bedroc`](api.md), [`calculate_rie`](api.md), [`calculate_enrichment_factor`](api.md), [`calculate_log_auc`](api.md), [`calculate_optimal_mcc`](api.md); all of them at once with [`evaluate_all_enrichment_metrics`](api.md).

**Where the thresholds come from.** They are conventional defaults for screening benchmarks, not universal limits. A virtual screen with AUC 0.65 can still be useful for a given target, and one with AUC 0.75 can be misleading. Adjust the thresholds to your use case and report them.

### Confidence intervals

- **ROC-AUC:** DeLong interval on the logit scale. In the binormal validation it covered the true AUC more often than the percentile bootstrap ([Validation](validation.md)).
- **Other metrics:** stratified non-parametric bootstrap (1,000 resamples by default, 95 % level, fixed seed), keeping the number of actives and decoys in each resample.

[`bootstrap_enrichment_ci`](api.md) returns the point estimate and the interval of each metric. With few actives the intervals are wide, and that is information: an enrichment factor from 10 actives is not precise.

## Redocking (pose accuracy)

The RMSD between each docked pose and the crystallographic ligand is computed on heavy atoms, with **symmetry correction** by graph isomorphism (equivalent atoms are not penalised, for example the two oxygens of a carboxylate) and **without superposition**, because redocking is done in the receptor frame. Functions: [`calculate_symmetry_corrected_rmsd`](api.md), [`calculate_heavy_atom_rmsd`](api.md).

[`evaluate_redocking_success`](api.md) judges the result:

| Input | Status |
|---|---|
| Ranked poses of one complex (`ranked_poses=True`) | `PASS` if the top-ranked pose is within 2 Å; `WARNING` if a lower-ranked pose is within 2 Å (a scoring failure) or the top-ranked pose is within 3 Å; `FAIL` otherwise |
| One RMSD per complex over a set (`ranked_poses=False`) | `PASS` if at least 70 % are within 2 Å; `WARNING` if at least 50 %; `FAIL` otherwise |

The 2 Å and 3 Å limits are arguments (`threshold_pass`, `threshold_warn`). The difference between "no pose reproduces the crystal pose" and "a pose does, but it is not ranked first" matters: the first is a sampling problem, the second a scoring problem.

## Decoy bias

If actives and decoys differ in simple properties such as molecular weight or logP, a model can separate them without capturing binding, and enrichment is inflated. When property columns are provided, DockCert compares the distributions of actives and decoys for each property (Kolmogorov–Smirnov and Wasserstein statistics) and counts the properties with a significant difference (\(p < 0.001\)).

| Properties with a significant difference | Risk | Status |
|---|---|---|
| 0 | `LOW` | `PASS` |
| 1 or 2 | `MODERATE` | `WARNING` |
| 3 or more | `HIGH` | `WARNING` |

The bias check raises a `WARNING` but never a `FAIL`, since the metrics themselves are not wrong; they may just not mean what they seem to. Function: [`evaluate_decoy_bias`](api.md).
