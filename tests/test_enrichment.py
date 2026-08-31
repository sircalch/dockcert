"""
Tests for virtual screening enrichment metrics (ROC-AUC, PR-AUC, BEDROC, RIE, EF).
"""

import numpy as np
import pytest
from dockcert.core.enrichment import (
    calculate_roc_auc,
    calculate_pr_auc,
    calculate_bedroc,
    calculate_rie,
    calculate_enrichment_factor,
    calculate_log_auc,
    calculate_optimal_mcc,
    evaluate_all_enrichment_metrics
)


def test_perfect_ranking():
    # 10 actives, 90 decoys, perfectly ranked first
    labels = np.array([1]*10 + [0]*90)
    scores = np.array(list(range(100)))  # 0 to 9 are best (lower is better)
    
    roc_auc = calculate_roc_auc(labels, scores, lower_is_better=True)
    assert np.isclose(roc_auc, 1.0)
    
    pr_auc = calculate_pr_auc(labels, scores, lower_is_better=True)
    assert np.isclose(pr_auc, 1.0)
    
    bedroc = calculate_bedroc(labels, scores, alpha=20.0, lower_is_better=True)
    assert np.isclose(bedroc, 1.0, atol=1e-3)
    
    ef1, ef1_max = calculate_enrichment_factor(labels, scores, fraction=0.01, lower_is_better=True)
    assert np.isclose(ef1, 10.0)  # 1/1 / (10/100) = 10.0
    
    mcc, thresh = calculate_optimal_mcc(labels, scores, lower_is_better=True)
    assert np.isclose(mcc, 1.0)


def test_random_ranking():
    rng = np.random.default_rng(42)
    n = 2000
    labels = rng.binomial(1, 0.05, size=n)
    scores = rng.normal(0, 1, size=n)
    
    roc_auc = calculate_roc_auc(labels, scores)
    assert 0.45 <= roc_auc <= 0.55
    
    bedroc = calculate_bedroc(labels, scores, alpha=20.0)
    assert 0.0 <= bedroc <= 0.15


def test_bedroc_bounds():
    labels = np.array([0]*90 + [1]*10)  # Actives at the very worst ranks
    scores = np.array(list(range(100)))
    
    bedroc = calculate_bedroc(labels, scores, alpha=20.0, lower_is_better=True)
    assert np.isclose(bedroc, 0.0, atol=1e-3)


def test_evaluate_all_metrics():
    labels = np.array([1]*5 + [0]*95)
    scores = np.linspace(-10, 0, 100)
    
    metrics_dict = evaluate_all_enrichment_metrics(labels, scores, lower_is_better=True)
    assert metrics_dict["roc_auc"] > 0.90
    assert metrics_dict["bedroc_20"] > 0.80
    assert metrics_dict["ef_1pct"] > 0.0
