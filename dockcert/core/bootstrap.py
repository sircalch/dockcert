"""
Stratified bootstrap confidence interval estimation for virtual screening metrics.
"""

from typing import Dict, Any, Tuple, Optional
import numpy as np
from dockcert.core.enrichment import (
    calculate_roc_auc,
    calculate_pr_auc,
    calculate_bedroc,
    calculate_enrichment_factor
)


def bootstrap_enrichment_ci(
    labels: np.ndarray,
    scores: np.ndarray,
    n_resamples: int = 1000,
    confidence_level: float = 0.95,
    lower_is_better: bool = True,
    random_state: Optional[int] = 42
) -> Dict[str, Tuple[float, float, float]]:
    """
    Computes stratified non-parametric bootstrap 95% Confidence Intervals for all key metrics.

    Parameters
    ----------
    labels : np.ndarray
        Binary label array (1 = Active, 0 = Decoy).
    scores : np.ndarray
        Docking scores.
    n_resamples : int, default 1000
        Number of bootstrap iterations.
    confidence_level : float, default 0.95
        Confidence level.
    lower_is_better : bool
        Sorting direction.
    random_state : int, optional
        Seed for reproducibility.

    Returns
    -------
    ci_dict : dict
        Mapping metric_name -> (point_estimate, ci_lower, ci_upper).
    """
    y_true = np.asarray(labels, dtype=int)
    y_scores = np.asarray(scores, dtype=float)
    
    active_idx = np.where(y_true == 1)[0]
    decoy_idx = np.where(y_true == 0)[0]
    
    n_act = len(active_idx)
    n_dec = len(decoy_idx)
    
    if n_act == 0 or n_dec == 0:
        return {}
        
    rng = np.random.default_rng(random_state)
    
    # Point estimates
    pe_roc = calculate_roc_auc(y_true, y_scores, lower_is_better)
    pe_pr = calculate_pr_auc(y_true, y_scores, lower_is_better)
    pe_bedroc20 = calculate_bedroc(y_true, y_scores, alpha=20.0, lower_is_better=lower_is_better)
    pe_bedroc80 = calculate_bedroc(y_true, y_scores, alpha=80.5, lower_is_better=lower_is_better)
    pe_ef1, _ = calculate_enrichment_factor(y_true, y_scores, fraction=0.01, lower_is_better=lower_is_better)
    pe_ef5, _ = calculate_enrichment_factor(y_true, y_scores, fraction=0.05, lower_is_better=lower_is_better)
    pe_ef10, _ = calculate_enrichment_factor(y_true, y_scores, fraction=0.10, lower_is_better=lower_is_better)
    
    boot_roc = np.empty(n_resamples, dtype=float)
    boot_pr = np.empty(n_resamples, dtype=float)
    boot_bedroc20 = np.empty(n_resamples, dtype=float)
    boot_bedroc80 = np.empty(n_resamples, dtype=float)
    boot_ef1 = np.empty(n_resamples, dtype=float)
    boot_ef5 = np.empty(n_resamples, dtype=float)
    boot_ef10 = np.empty(n_resamples, dtype=float)
    
    for b in range(n_resamples):
        sample_act = rng.choice(active_idx, size=n_act, replace=True)
        sample_dec = rng.choice(decoy_idx, size=n_dec, replace=True)
        
        sample_idx = np.concatenate([sample_act, sample_dec])
        b_labels = y_true[sample_idx]
        b_scores = y_scores[sample_idx]
        
        boot_roc[b] = calculate_roc_auc(b_labels, b_scores, lower_is_better)
        boot_pr[b] = calculate_pr_auc(b_labels, b_scores, lower_is_better)
        boot_bedroc20[b] = calculate_bedroc(b_labels, b_scores, alpha=20.0, lower_is_better=lower_is_better)
        boot_bedroc80[b] = calculate_bedroc(b_labels, b_scores, alpha=80.5, lower_is_better=lower_is_better)
        ef1, _ = calculate_enrichment_factor(b_labels, b_scores, fraction=0.01, lower_is_better=lower_is_better)
        ef5, _ = calculate_enrichment_factor(b_labels, b_scores, fraction=0.05, lower_is_better=lower_is_better)
        ef10, _ = calculate_enrichment_factor(b_labels, b_scores, fraction=0.10, lower_is_better=lower_is_better)
        boot_ef1[b] = ef1
        boot_ef5[b] = ef5
        boot_ef10[b] = ef10
        
    alpha_ci = (1.0 - confidence_level) / 2.0
    
    def get_ci(point_est, boot_arr):
        low = float(np.percentile(boot_arr, 100.0 * alpha_ci))
        high = float(np.percentile(boot_arr, 100.0 * (1.0 - alpha_ci)))
        return float(point_est), low, high

    return {
        "roc_auc": get_ci(pe_roc, boot_roc),
        "pr_auc": get_ci(pe_pr, boot_pr),
        "bedroc_20": get_ci(pe_bedroc20, boot_bedroc20),
        "bedroc_80": get_ci(pe_bedroc80, boot_bedroc80),
        "ef_1pct": get_ci(pe_ef1, boot_ef1),
        "ef_5pct": get_ci(pe_ef5, boot_ef5),
        "ef_10pct": get_ci(pe_ef10, boot_ef10)
    }
