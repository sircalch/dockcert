# Validation

The `validation/` folder of the repository contains the scripts and the result files behind these benchmarks, so they can be repeated.

| Question | Data | Result | Files |
|---|---|---|---|
| Does the RMSD match an independent implementation? | 85 Astex Diverse Set complexes, redocked with AutoDock Vina (9 poses each) | Top-1 success (RMSD \(\le 2\) Å): **61.2 %** (95 % CI 50.0 to 71.6 %), identical to RDKit `CalcRMS`; success in the best of 9 poses: 82.4 % | `astex_redock.py`, `analyze_astex.py`, `results/astex_success.csv` |
| Does symmetry correction matter? | Same 85 complexes | Without symmetry correction, top-1 success drops to 55.3 %; version 1.0.0 (a per-element Hungarian assignment) gave 20.0 % | `results/astex_success.csv` |
| Does the AUC interval reach its nominal level? | Binormal simulations with known AUC (20 to 100 actives) | Coverage of the 95 % interval: 0.89 to 0.95 for the DeLong logit interval, 0.87 to 0.95 for the percentile bootstrap | `auc_ci_coverage.py`, `bootstrap_coverage.py`, `results/auc_ci_coverage.csv`, `results/bootstrap_coverage.csv` |
| What does a real screen look like? | DUD-E targets AmpC (48 actives, 2,850 decoys) and InhA (43 actives, 2,300 decoys) with Vina scores | AUC 0.60 (95 % CI 0.52 to 0.67) for AmpC and 0.63 (0.56 to 0.71) for InhA | `dude_screen.py`, `analyze_dude.py`, `results/dude_metrics.csv` |

## Reading the results with care

- **The 1.0.0 to 1.1.0 change in RMSD** is the clearest finding: the earlier assignment method undercounted successes by a large factor on real data, which is why the validation compares against an independent implementation.
- **The DUD-E results are modest.** With default Vina scoring, enrichment on these two targets is only slightly better than random. DockCert reports that as it is; the thresholds in [How the checks work](quality-control.md) would flag these screens as marginal or failing, which is the intended behaviour.
- **Confidence intervals are not perfect.** With few actives, the coverage of both AUC intervals falls to about 0.87 to 0.91 in some settings.

The simulation settings, seeds and software versions are in the scripts and result files of `validation/`.
