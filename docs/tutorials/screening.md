# Screening enrichment from a table

## The input table

A CSV (or TSV) with one row per compound, a docking score and an active/decoy label. The score and label columns are detected from their names; otherwise name them with `--score-col` and `--label-col`. Other numeric columns, such as molecular weight or logP, are used for the decoy-bias check.

The repository includes `benchmark_screening.csv`, a synthetic set of 1,550 compounds (50 actives, 1,500 decoys):

```text
Compound_ID,Docking_Score_kcal_mol,Activity_Status,Molecular_Weight,LogP
DECOY_1470,-7.600870373071028,decoy,441.92676327844185,3.7169610482990274
DECOY_1392,-6.233458460592504,decoy,344.95163035645083,2.5226916288742673
```

By default **lower scores are better**, as for AutoDock Vina binding energies. For scoring functions where higher is better, add `--higher-is-better` (Python: `lower_is_better=False`).

## Command line

```bash
dockcert assess -i benchmark_screening.csv -o screening_report/
```

```text
 [RESULT] Overall Docking Certification: PASS
 * ROC-AUC           :  0.973 | Status: PASS
 * BEDROC (alpha=20.0):  0.784 | Status: PASS
 * EF 1%             : 29.062 | Status: PASS
 * EF 5%             : 14.705 | Status: PASS
 * Optimal MCC       :  0.610 | Status: PASS
```

Open `screening_report/report.html` for the confidence intervals, the ROC and enrichment curves, and the diagnostic message of each check.

## Python API

```python
from dockcert import assess_docking_quality
from dockcert.parsers.generic_csv import load_docking_csv

labels, scores, properties, meta = load_docking_csv("benchmark_screening.csv")

report = assess_docking_quality(
    labels=labels,
    scores=scores,
    active_properties={k: v[labels == 1] for k, v in properties.items()},
    decoy_properties={k: v[labels == 0] for k, v in properties.items()},
    lower_is_better=True,
)
print(report.overall_status)
for item in report.enrichment_metrics.values():
    print(item.name, round(item.value, 3), item.status)
```

Passing `active_properties` and `decoy_properties` turns on the [decoy-bias check](../quality-control.md#decoy-bias). To also judge pose accuracy, see the [worked example](worked-example.md).

## Reading the result

A `PASS` for every metric on a synthetic set like this one is expected. On real targets the numbers are usually more modest; see [Validation](../validation.md) for DUD-E results where AUC is around 0.6.
