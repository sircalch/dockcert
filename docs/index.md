# DockCert

**Statistical validation and reproducibility assessment for molecular docking and virtual screening.**

A docking study usually reports one number: an enrichment factor, an AUC, or the fraction of poses within 2 Å. DockCert makes the evidence behind that number explicit. It reads a screening table, optionally a crystallographic ligand and docked poses, and checks:

- **Virtual screening performance:** ROC-AUC, precision–recall AUC, BEDROC, RIE, enrichment factors, log AUC and the optimal Matthews correlation coefficient, with confidence intervals and a pass/warning/fail judgement against explicit thresholds.
- **Pose accuracy (redocking):** heavy-atom RMSD to the crystallographic ligand with symmetry correction, and whether the top-ranked pose or only a lower-ranked one reproduces the binding mode.
- **Decoy bias:** whether actives and decoys differ in physicochemical properties such as molecular weight and logP, which can inflate enrichment artificially.

Every decision is recorded, and the results are written as a self-contained HTML report, figures, LaTeX and CSV tables and a draft methods paragraph.

```bash
pip install dockcert
dockcert demo -o my_docking_validation/
```

Open `my_docking_validation/report.html` in a browser to see a report built from a synthetic DUD-E-like dataset.

## Where to go next

- [Installation](installation.md)
- [Screening enrichment from a table](tutorials/screening.md)
- [Worked example: redocking with AutoDock Vina](tutorials/worked-example.md)
- [How the checks work](quality-control.md): metrics, thresholds and what each status means
- [Validation](validation.md): benchmarks on 85 redocking cases and two DUD-E targets
- [API reference](api.md)

## Scope

DockCert checks the **evaluation** of docking results. It does not run docking, does not score ligands and does not replace expert judgement: a `PASS` means the checks implemented here found no problem with the evidence, not that the docking protocol is right for your target.

## Getting help

Questions and bug reports are welcome in the [issue tracker](https://github.com/sircalch/dockcert/issues). See [Support and contributing](support.md).
