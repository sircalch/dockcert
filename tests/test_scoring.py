"""
Tests for scoring decision matrix, report generation, and CLI demo in DockCert.
"""

import os
import tempfile
import numpy as np
import pytest
from dockcert.core.scoring import assess_docking_quality
from dockcert.reporters.plot_generator import generate_docking_figures
from dockcert.reporters.manuscript_prep import generate_docking_manuscript_assets
from dockcert.reporters.html_report import generate_docking_html_report
from dockcert.cli import run_demo


def test_full_docking_validation_pipeline():
    rng = np.random.default_rng(42)
    n_act = 30
    n_dec = 300
    
    act_scores = rng.normal(-9.5, 0.8, n_act)
    dec_scores = rng.normal(-6.5, 1.0, n_dec)
    
    labels = np.concatenate([np.ones(n_act), np.zeros(n_dec)])
    scores = np.concatenate([act_scores, dec_scores])
    
    report = assess_docking_quality(
        labels=labels,
        scores=scores,
        rmsd_values=[1.35, 1.80, 2.50],
        lower_is_better=True
    )
    
    assert report.overall_status == "PASS"
    assert "roc_auc" in report.enrichment_metrics
    assert report.enrichment_metrics["roc_auc"].value > 0.80
    assert report.redocking_result["status"] == "PASS"
    
    with tempfile.TemporaryDirectory() as tmpdir:
        # Plot generation
        plots = generate_docking_figures(labels, scores, report, tmpdir, formats=["png", "svg"])
        assert len(plots) > 0
        for p in plots:
            assert os.path.exists(p)
            
        # Manuscript assets
        assets = generate_docking_manuscript_assets(report, tmpdir)
        assert os.path.exists(assets["summary_csv"])
        assert os.path.exists(assets["summary_tex"])
        assert os.path.exists(assets["methods_text"])
        assert os.path.exists(assets["citation_bib"])
        
        # HTML report
        html_p = os.path.join(tmpdir, "report.html")
        generate_docking_html_report(report, html_p, methods_text="Sample", citation_bib="@software{}")
        assert os.path.exists(html_p)
        assert os.path.getsize(html_p) > 500


def test_cli_demo_execution():
    with tempfile.TemporaryDirectory() as tmpdir:
        run_demo(output_dir=tmpdir)
        assert os.path.exists(os.path.join(tmpdir, "report.html"))
        assert os.path.exists(os.path.join(tmpdir, "dockcert_summary_table.csv"))
        assert os.path.exists(os.path.join(tmpdir, "citation.bib"))
