"""
Quickstart API tutorial for DockCert.
"""

import os
import numpy as np
from dockcert import assess_docking_quality
from dockcert.reporters import (
    generate_docking_figures,
    generate_docking_manuscript_assets,
    generate_docking_html_report
)


def main():
    print("Running DockCert Python API quickstart example...")
    output_dir = "quickstart_dockcert_output"
    os.makedirs(output_dir, exist_ok=True)
    
    # 1. Simulate 40 actives and 800 decoys
    rng = np.random.default_rng(42)
    n_act = 40
    n_dec = 800
    
    act_scores = rng.normal(-9.5, 0.9, n_act)
    dec_scores = rng.normal(-6.6, 1.0, n_dec)
    
    labels = np.concatenate([np.ones(n_act, dtype=int), np.zeros(n_dec, dtype=int)])
    scores = np.concatenate([act_scores, dec_scores])
    
    # Redocking RMSD values for top 5 poses
    redocking_rmsds = [1.38, 1.75, 2.10, 2.80, 3.40]
    
    # 2. Assess Quality
    report = assess_docking_quality(
        labels=labels,
        scores=scores,
        rmsd_values=redocking_rmsds,
        lower_is_better=True
    )
    
    print(f"\nOverall Docking Certification: {report.overall_status}")
    print(f"Validation Score: {report.validation_score}")
    for k, item in report.enrichment_metrics.items():
        print(f" - {item.name:18s}: {item.value:.3f} | Status: {item.status}")
    if report.redocking_result:
        print(f" - Redocking Best RMSD: {report.redocking_result['min_rmsd']:.2f} A | Status: {report.redocking_result['status']}")
        
    # 3. Export all assets
    generate_docking_figures(labels, scores, report, output_dir, lower_is_better=True, rmsd_values=redocking_rmsds)
    assets = generate_docking_manuscript_assets(report, output_dir)
    
    with open(assets["methods_text"], "r", encoding="utf-8") as f:
        methods_txt = f.read()
    with open(assets["citation_bib"], "r", encoding="utf-8") as f:
        bib_txt = f.read()
        
    html_p = os.path.join(output_dir, "report.html")
    generate_docking_html_report(report, html_p, methods_text=methods_txt, citation_bib=bib_txt)
    
    print(f"\nCompleted! Check out: {os.path.abspath(html_p)}")


if __name__ == "__main__":
    main()
