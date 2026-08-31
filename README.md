# DockCert

[![CI](https://github.com/amonreal/dockcert/actions/workflows/test.yml/badge.svg)](https://github.com/amonreal/dockcert/actions)
[![PyPI version](https://img.shields.io/pypi/v/dockcert.svg?color=blue)](https://pypi.org/project/dockcert/)
[![Python versions](https://img.shields.io/pypi/pyversions/dockcert.svg)](https://pypi.org/project/dockcert/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1234568.svg)](https://doi.org/10.5281/zenodo.1234568)

> **Automated Statistical Validation, Early Enrichment Metrics, and Reproducibility Assessment for Molecular Docking and Virtual Screening Studies.**

---

## Overview

**DockCert** is an open-source scientific toolkit engineered to evaluate whether a molecular docking or virtual screening protocol is rigorously validated before publication.

Instead of computing enrichment factors or pose RMSD across fragmented scripts, `dockcert` analyzes your screening results and docked structures in a single run:

- 🎯 **Redocking & Cross-Docking Pose Accuracy**: Heavy-atom & symmetry-corrected RMSD (Hungarian algorithm).
- 📈 **Early Enrichment Quantification**:
  - **BEDROC** ($\alpha=20.0, 80.5, 160.9$)
  - **RIE** (Robust Initial Enhancement)
  - **EF1%**, **EF5%**, **EF10%** (Enrichment Factors)
  - **ROC-AUC** & **PR-AUC** (Precision-Recall AUC)
  - **logAUC** & **Optimal MCC** (Matthews Correlation Coefficient)
- 🎲 **95% Stratified Bootstrap Confidence Intervals** for every metric.
- 🔬 **Decoy Bias & Artificial Enrichment Audit**: Kolmogorov-Smirnov & Wasserstein property tests (MW, LogP, HBD, HBA).
- 🚦 **Quality Certification Badges (`PASS` / `WARNING` / `FAIL`)**.
- 📑 **Publication Deliverables**: Interactive self-contained `report.html`, vector plots (SVG/PDF/PNG 300 DPI), LaTeX summary tables (`.tex`), and a draft **Methods & Supporting Information** paragraph with automated **BibTeX citations**.
- 🔌 **Engine Agnostic**: Works seamlessly with AutoDock Vina, Smina, GNINA, Schrödinger Glide, CCDC GOLD, DOCK, rDock, and generic CSV tables.

```
  Docking Results (.csv, .sdf, .pdbqt, logs)
                     │
                     ▼
  ┌───────────────────────────────────────────────────────────┐
  │                         DockCert                          │
  │  ├── Redocking RMSD (Symmetry corrected)                  │
  │  ├── Early Enrichment (BEDROC, RIE, EF1%, ROC-AUC)        │
  │  ├── Stratified 95% Bootstrap Confidence Intervals        │
  │  └── Decoy Bias Audit (MW, LogP KS-tests)                 │
  └───────────────────────────────────────────────────────────┘
                     │
                     ▼
  ┌───────────────────────────────────────────────────────────┐
  │                   Publication Deliverables                │
  │  ├── report.html (Interactive Dashboard & Badges)         │
  │  ├── dockcert_validation_overview.pdf/svg/png             │
  │  ├── dockcert_summary_table.tex / .csv                    │
  │  ├── methods_snippet.txt (Ready for Manuscript)           │
  │  └── citation.bib (BibTeX Reference)                      │
  └───────────────────────────────────────────────────────────┘
```

---

## Installation

### From PyPI
```bash
pip install dockcert
```

### From Source
```bash
git clone https://github.com/amonreal/dockcert.git
cd dockcert
pip install -e .[dev]
```

---

## Quickstart (CLI)

### 1. Run Demonstration Mode (Instant Benchmark Dataset)
```bash
dockcert demo -o my_docking_validation/
```
Open `my_docking_validation/report.html` in your browser!

### 2. Assess Virtual Screening CSV
```bash
dockcert assess -i screening_results.csv --score-col docking_score --label-col is_active -o validation_report/
```

### 3. Assess with Reference Ligand Pose RMSD
```bash
dockcert assess -i screening_results.csv --ref-ligand crystal_ligand.sdf --docked-pose docked_pose.sdf -o full_report/
```

---

## Python API Usage

```python
import numpy as np
from dockcert import assess_docking_quality
from dockcert.reporters import generate_docking_figures, generate_docking_manuscript_assets, generate_docking_html_report

# Screening scores (e.g. 50 actives, 1000 decoys)
labels = np.array([1]*50 + [0]*1000)
scores = np.concatenate([np.random.normal(-9.2, 0.8, 50), np.random.normal(-6.5, 1.1, 1000)])

# Assess quality
report = assess_docking_quality(
    labels=labels,
    scores=scores,
    rmsd_values=[1.42, 1.85, 2.30],
    lower_is_better=True
)

print(f"Overall Validation Status: {report.overall_status}")
print(f"ROC-AUC: {report.enrichment_metrics['roc_auc'].value:.3f}")
print(f"BEDROC (alpha=20.0): {report.enrichment_metrics['bedroc_20'].value:.3f}")
print(f"EF 1%: {report.enrichment_metrics['ef_1pct'].value:.1f}x")

# Export publication deliverables
generate_docking_figures(labels, scores, report, "output_dir/")
generate_docking_manuscript_assets(report, "output_dir/")
generate_docking_html_report(report, "output_dir/report.html")
```

---

## Citation

If you use DockCert to evaluate molecular docking validation or virtual screening enrichment, please cite:

```bibtex
@software{monreal2026dockcert,
  author = {Monreal-Hern{\'a}ndez, Andre},
  title = {{DockCert: An Open-Source Toolkit for Statistical Validation, Enrichment Metrics, and Reproducibility Assessment of Molecular Docking Studies}},
  year = {2026},
  version = {1.0.0},
  publisher = {Zenodo},
  url = {https://github.com/amonreal/dockcert}
}
```

---

## License

This project is licensed under the **MIT License** - see the [LICENSE](LICENSE) file for details.
