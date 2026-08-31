"""
Manuscript Methods text generator, LaTeX summary tables, and BibTeX citations for DockCert.
"""

from typing import Dict, Any
import os
import pandas as pd
from dockcert.core.scoring import DockingValidationReport


def generate_docking_manuscript_assets(
    report: DockingValidationReport,
    output_dir: str
) -> Dict[str, str]:
    """
    Generates all textual and tabular assets needed for publication manuscripts.

    Parameters
    ----------
    report : DockingValidationReport
        Docking validation report.
    output_dir : str
        Target output directory.

    Returns
    -------
    paths : dict
        Mapping asset_key -> file path.
    """
    os.makedirs(output_dir, exist_ok=True)
    generated = {}
    
    # 1. Summary DataFrame
    rows = []
    
    if report.redocking_result is not None:
        rr = report.redocking_result
        rows.append({
            "Validation Metric": "Redocking Pose RMSD (Best)",
            "Value": f"{rr['min_rmsd']:.2f} A",
            "95% Bootstrap CI": "N/A",
            "Pass Threshold": "<= 2.0 A",
            "Status": rr["status"],
            "Diagnostic Finding": rr["recommendation"]
        })
        
    for key, item in report.enrichment_metrics.items():
        ci_str = f"[{item.ci_lower_95:.3f}, {item.ci_upper_95:.3f}]" if item.ci_lower_95 is not None else "N/A"
        th_str = f">= {item.threshold_pass:.2f}"
        rows.append({
            "Validation Metric": item.name,
            "Value": f"{item.value:.3f}",
            "95% Bootstrap CI": ci_str,
            "Pass Threshold": th_str,
            "Status": item.status,
            "Diagnostic Finding": item.message
        })
        
    if report.decoy_bias_result is not None:
        db = report.decoy_bias_result
        rows.append({
            "Validation Metric": "Decoy Property Bias",
            "Value": f"Risk: {db['risk_level']}",
            "95% Bootstrap CI": "N/A",
            "Pass Threshold": "LOW Risk",
            "Status": db["status"],
            "Diagnostic Finding": db["recommendation"]
        })
        
    df_summary = pd.DataFrame(rows)
    
    # CSV Table
    csv_path = os.path.join(output_dir, "dockcert_summary_table.csv")
    df_summary.to_csv(csv_path, index=False)
    generated["summary_csv"] = csv_path
    
    # LaTeX Table
    tex_table_path = os.path.join(output_dir, "dockcert_summary_table.tex")
    tex_table = df_summary.to_latex(index=False, escape=False)
    with open(tex_table_path, "w", encoding="utf-8") as f:
        f.write("% DockCert Molecular Docking & Virtual Screening Statistical Validation Table\n")
        f.write(tex_table)
    generated["summary_tex"] = tex_table_path

    # 2. Methods Text Snippet
    methods_path = os.path.join(output_dir, "methods_snippet.txt")
    
    redock_text = ""
    if report.redocking_result is not None:
        best_rmsd = report.redocking_result["min_rmsd"]
        redock_text = (
            f"The docking protocol was structurally validated by redocking the co-crystallized reference ligand, "
            f"achieving a heavy-atom Root-Mean-Square Deviation (RMSD) of {best_rmsd:.2f} \\AA relative to the experimental pose (status: {report.redocking_result['status']}). "
        )
        
    enrich_text = ""
    if report.enrichment_metrics:
        roc_val = report.enrichment_metrics.get("roc_auc")
        bedroc_val = report.enrichment_metrics.get("bedroc_20")
        ef1_val = report.enrichment_metrics.get("ef_1pct")
        
        roc_str = f"ROC-AUC = {roc_val.value:.3f} [{roc_val.ci_lower_95:.2f}, {roc_val.ci_upper_95:.2f}]" if roc_val and roc_val.ci_lower_95 else (f"ROC-AUC = {roc_val.value:.3f}" if roc_val else "")
        bedroc_str = f"BEDROC (\\alpha = 20.0) = {bedroc_val.value:.3f}" if bedroc_val else ""
        ef1_str = f"EF_{{1\\%}} = {ef1_val.value:.1f}" if ef1_val else ""
        
        metrics_combined = ", ".join([s for s in [roc_str, bedroc_str, ef1_str] if s])
        
        bias_str = ""
        if report.decoy_bias_result is not None:
            bias_str = f" Physicochemical decoy matching was verified with {report.decoy_bias_result['risk_level'].lower()} risk of artificial bias."
            
        enrich_text = (
            f"Virtual screening discrimination and early enrichment were statistically assessed using DockCert v1.0.0 (Monreal-Hernández, 2026), "
            f"yielding {metrics_combined} with 95% stratified bootstrap confidence intervals across {report.dataset_summary.get('n_actives', 0)} active and {report.dataset_summary.get('n_decoys', 0)} decoy molecules.{bias_str}"
        )
        
    full_methods = (
        f"{redock_text}{enrich_text} Overall docking validation status was certified as {report.overall_status}."
    )
    
    with open(methods_path, "w", encoding="utf-8") as f:
        f.write(full_methods + "\n")
    generated["methods_text"] = methods_path

    # 3. BibTeX Citation File
    bib_path = os.path.join(output_dir, "citation.bib")
    bib_content = """@software{monreal2026dockcert,
  author = {Monreal-Hern\\'andez, Andre},
  title = {{DockCert: An Open-Source Toolkit for Statistical Validation, Enrichment Metrics, and Reproducibility Assessment of Molecular Docking Studies}},
  year = {2026},
  version = {1.0.0},
  publisher = {Zenodo},
  url = {https://github.com/amonreal/dockcert}
}
"""
    with open(bib_path, "w", encoding="utf-8") as f:
        f.write(bib_content)
    generated["citation_bib"] = bib_path

    return generated
