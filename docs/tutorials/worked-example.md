# Worked example: redocking with AutoDock Vina

This example uses real results: AutoDock Vina 1.2.7 redocking of three crystal complexes, whose files are bundled with the repository in `tests/data/redock/`:

| Complex | System |
|---|---|
| 3PTB | trypsin with benzamidine |
| 1STP | streptavidin with biotin |
| 1HSG | HIV-1 protease with indinavir |

For each, the folder has the crystallographic ligand (`<PDB>_ref.sdf`), the nine ranked poses from Vina (`<PDB>_out.pdbqt`) and the Vina log. Run the example from the repository root:

```bash
python examples/redock_to_report.py
```

## 1. Does the top-ranked pose reproduce the crystal pose?

```text
complex   top-1 (A)  best of 9 (A)  status
3PTB           0.38           0.38  PASS   (trypsin-benzamidine)
1STP           0.60           0.60  PASS   (streptavidin-biotin)
1HSG           4.39           4.11  FAIL   (HIV-1 protease-indinavir)
```

DockCert computes the heavy-atom RMSD of every pose to the crystallographic ligand with **symmetry correction** (a graph-isomorphism match, so equivalent atoms are not penalised) and without superposing the poses, because redocking happens in the receptor frame.

- **3PTB and 1STP** reproduce the crystal pose with the top-ranked pose (below 2 Å).
- **1HSG** fails: no pose of the nine is within 2 Å (best 4.11 Å). Indinavir is a large, flexible ligand with 45 heavy atoms, a typical case where a docking protocol does not reproduce the binding mode. DockCert reports this as a failure of the protocol for this target rather than hiding it in an average.

The same checks from the command line:

```bash
dockcert assess -i benchmark_screening.csv \
    --ref-ligand tests/data/redock/1HSG_ref.sdf \
    --docked-pose tests/data/redock/1HSG_out.pdbqt \
    --vina-log tests/data/redock/1HSG_vina.log -o report_1HSG/
```

```text
 [RESULT] Overall Docking Certification: FAIL
 * Redocking RMSD     : top-1 4.39 A, best of 9 4.11 A | Status: FAIL
```

## 2. Screening enrichment

The script also assesses the bundled screening table together with the redocking result of 3PTB:

```text
Overall status: PASS
 - ROC-AUC               0.973  [0.954, 0.984]     PASS
 - BEDROC (alpha=20.0)    0.784  [0.701, 0.855]     PASS
 - EF 1%                29.062  [25.188, 31.000]   PASS
 - EF 5%                14.705  [12.321, 16.692]   PASS
 - Optimal MCC           0.610                     PASS
```

The ROC-AUC interval is a DeLong interval on the logit scale and the others are bootstrap intervals; see [How the checks work](../quality-control.md). The overall status is the worst of all checks, so a protocol that fails redocking cannot be certified by enrichment alone.

## 3. The report

The script writes the files described in [Outputs](../outputs.md). `report.html` shows each metric with its interval and threshold, the redocking RMSDs and the decoy-bias audit. `methods_snippet.txt` contains a draft methods paragraph; edit it so that it describes your study.
