# Outputs

`dockcert assess` and `dockcert demo` write these files to the output directory:

| File | Content |
|---|---|
| `report.html` | Self-contained report with every check, its status, thresholds and diagnostics |
| `dockcert_validation_overview.svg`, `.pdf`, `.png` | ROC curve, enrichment and metric summary (when screening data is given) |
| `dockcert_redocking_rmsd.svg`, `.pdf`, `.png` | RMSD of the docked poses (when redocking data is given) |
| `dockcert_summary_table.csv`, `.tex` | Summary table, ready for a spreadsheet or a LaTeX manuscript |
| `methods_snippet.txt` | Draft methods paragraph describing the checks and their results |
| `citation.bib` | BibTeX entry for DockCert |

The methods paragraph is a starting point: edit it so that it describes your study accurately.
