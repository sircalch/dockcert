"""
Virtual Screening Enrichment Metrics: ROC-AUC, PR-AUC, BEDROC, RIE, EF, logAUC, and MCC.
"""

from typing import Dict, Any, Optional, Tuple, List
import numpy as np
from scipy import integrate
from sklearn import metrics


def _trapezoid_integral(y: np.ndarray, x: np.ndarray) -> float:
    """Helper to compute trapezoid integral across numpy 1.x, 2.x and scipy versions."""
    if hasattr(np, "trapezoid"):
        return float(np.trapezoid(y, x))
    elif hasattr(integrate, "trapezoid"):
        return float(integrate.trapezoid(y, x))
    elif hasattr(integrate, "trapz"):
        return float(integrate.trapz(y, x))
    elif hasattr(np, "trapz"):
        return float(np.trapz(y, x))
    else:
        # Fallback implementation
        return float(np.sum(0.5 * (y[:-1] + y[1:]) * np.diff(x)))


def calculate_roc_auc(
    labels: np.ndarray,
    scores: np.ndarray,
    lower_is_better: bool = True
) -> float:
    """
    Computes the Receiver Operating Characteristic Area Under the Curve (ROC-AUC).

    Parameters
    ----------
    labels : np.ndarray
        Binary array (1 = Active, 0 = Decoy/Inactive).
    scores : np.ndarray
        Docking scores or binding affinities.
    lower_is_better : bool, default True
        If True (standard for docking energies in kcal/mol), lower scores indicate higher affinity.

    Returns
    -------
    auc : float
        ROC-AUC value in [0, 1].
    """
    y_true = np.asarray(labels, dtype=int)
    y_scores = -np.asarray(scores, dtype=float) if lower_is_better else np.asarray(scores, dtype=float)
    
    if len(np.unique(y_true)) < 2:
        return 0.5
        
    auc = metrics.roc_auc_score(y_true, y_scores)
    return float(auc)


def calculate_pr_auc(
    labels: np.ndarray,
    scores: np.ndarray,
    lower_is_better: bool = True
) -> float:
    """
    Computes the Precision-Recall Area Under the Curve (PR-AUC / Average Precision).

    Parameters
    ----------
    labels : np.ndarray
        Binary array (1 = Active, 0 = Decoy).
    scores : np.ndarray
        Docking scores.
    lower_is_better : bool, default True
        If True, lower scores rank higher.

    Returns
    -------
    pr_auc : float
        PR-AUC in [0, 1].
    """
    y_true = np.asarray(labels, dtype=int)
    y_scores = -np.asarray(scores, dtype=float) if lower_is_better else np.asarray(scores, dtype=float)
    
    if len(np.unique(y_true)) < 2:
        return float(np.mean(y_true))
        
    return float(metrics.average_precision_score(y_true, y_scores))


def calculate_rie(
    labels: np.ndarray,
    scores: np.ndarray,
    alpha: float = 20.0,
    lower_is_better: bool = True
) -> float:
    """
    Computes the Robust Initial Enhancement (RIE) according to Truchon & Bayly (2007).

    Parameters
    ----------
    labels : np.ndarray
        Binary array (1 = Active, 0 = Decoy).
    scores : np.ndarray
        Docking scores.
    alpha : float, default 20.0
        Exponential weight factor (e.g. 20.0 weights the top 8% of the database).
    lower_is_better : bool, default True
        Sort direction.

    Returns
    -------
    rie : float
        Robust Initial Enhancement metric.
    """
    y_true = np.asarray(labels, dtype=int)
    y_scores = np.asarray(scores, dtype=float)
    
    order = np.argsort(y_scores) if lower_is_better else np.argsort(-y_scores)
    sorted_labels = y_true[order]
    
    n_total = len(sorted_labels)
    n_actives = int(np.sum(sorted_labels))
    
    if n_actives == 0 or n_actives == n_total or n_total < 2:
        return 1.0
        
    active_ranks = np.where(sorted_labels == 1)[0] + 1  # 1-indexed ranks
    
    sum_exp = np.sum(np.exp(-alpha * active_ranks / n_total))
    ra = (n_actives / n_total) * ((1.0 - np.exp(-alpha)) / (np.exp(alpha / n_total) - 1.0))
    
    if ra == 0:
        return 1.0
        
    rie = sum_exp / ra
    return float(rie)


def calculate_bedroc(
    labels: np.ndarray,
    scores: np.ndarray,
    alpha: float = 20.0,
    lower_is_better: bool = True
) -> float:
    """
    Computes the Boltzmann-Enhanced Discrimination of ROC (BEDROC) metric in [0, 1].

    Formula (Truchon & Bayly, J. Chem. Inf. Model. 2007, 47, 488-508):
        BEDROC = (S - S_min) / (S_max - S_min)
    where S = sum_{i=1}^n exp(-alpha * r_i / N)

    Parameters
    ----------
    labels : np.ndarray
        Binary array (1 = Active, 0 = Decoy).
    scores : np.ndarray
        Docking scores.
    alpha : float, default 20.0
        Exponential weight parameter (alpha=20.0 for top 8%, alpha=80.5 for top 2%).
    lower_is_better : bool, default True
        Sorting order.

    Returns
    -------
    bedroc : float
        Bounded early enrichment score in [0, 1].
    """
    y_true = np.asarray(labels, dtype=int)
    y_scores = np.asarray(scores, dtype=float)
    
    order = np.argsort(y_scores) if lower_is_better else np.argsort(-y_scores)
    sorted_labels = y_true[order]
    
    n_total = len(sorted_labels)
    n_actives = int(np.sum(sorted_labels))
    
    if n_actives == 0 or n_total < 2:
        return 0.0
    if n_actives == n_total:
        return 1.0
        
    active_ranks = np.where(sorted_labels == 1)[0] + 1  # 1-indexed ranks
    
    # Observed sum of exponential weights
    s = np.sum(np.exp(-alpha * active_ranks / n_total))
    
    # Theoretical maximum and minimum sums of exponential weights
    # S_max: ranks 1, 2, ..., n_actives
    # S_min: ranks N - n_actives + 1, ..., N
    s_max = np.sum(np.exp(-alpha * np.arange(1, n_actives + 1) / n_total))
    s_min = np.sum(np.exp(-alpha * np.arange(n_total - n_actives + 1, n_total + 1) / n_total))
    
    denom = s_max - s_min
    if denom == 0:
        return 1.0
        
    bedroc = (s - s_min) / denom
    return float(np.clip(bedroc, 0.0, 1.0))


def calculate_enrichment_factor(
    labels: np.ndarray,
    scores: np.ndarray,
    fraction: float = 0.01,
    lower_is_better: bool = True
) -> Tuple[float, float]:
    """
    Calculates the Enrichment Factor (EF_x%) and the theoretical maximum EF.

    Formula:
        EF = (Actives_top_x / N_top_x) / (Actives_total / N_total)

    Parameters
    ----------
    labels : np.ndarray
        Binary array (1 = Active, 0 = Decoy).
    scores : np.ndarray
        Docking scores.
    fraction : float, default 0.01
        Top database fraction to evaluate (e.g. 0.01 for EF1%, 0.05 for EF5%).
    lower_is_better : bool, default True
        Sorting order.

    Returns
    -------
    ef : float
        Observed Enrichment Factor.
    ef_max : float
        Theoretical maximum possible EF for the dataset at this fraction.
    """
    y_true = np.asarray(labels, dtype=int)
    y_scores = np.asarray(scores, dtype=float)
    
    n_total = len(y_true)
    n_actives = int(np.sum(y_true))
    
    if n_actives == 0 or n_total == 0:
        return 0.0, 0.0
        
    order = np.argsort(y_scores) if lower_is_better else np.argsort(-y_scores)
    sorted_labels = y_true[order]
    
    cutoff_n = max(1, int(np.round(fraction * n_total)))
    actives_in_top = np.sum(sorted_labels[:cutoff_n])
    
    # Enrichment Factor
    ef = (actives_in_top / cutoff_n) / (n_actives / n_total)
    
    # Maximum possible EF
    ef_max = min(1.0 / fraction, float(n_total) / float(n_actives))
    
    return float(ef), float(ef_max)


def calculate_log_auc(
    labels: np.ndarray,
    scores: np.ndarray,
    lam: float = 0.001,
    lower_is_better: bool = True
) -> float:
    """
    Computes semilogarithmic ROC area (logAUC) focusing on early False Positive Rates.

    Parameters
    ----------
    labels : np.ndarray
        Binary array.
    scores : np.ndarray
        Docking scores.
    lam : float, default 0.001
        Lower integration limit on FPR (e.g. 0.1% = 0.001).
    lower_is_better : bool
        Sorting order.

    Returns
    -------
    log_auc : float
        logAUC value in [0, 1].
    """
    y_true = np.asarray(labels, dtype=int)
    y_scores = -np.asarray(scores, dtype=float) if lower_is_better else np.asarray(scores, dtype=float)
    
    if len(np.unique(y_true)) < 2:
        return 0.5
        
    fpr, tpr, _ = metrics.roc_curve(y_true, y_scores)
    
    # Filter for fpr >= lam
    fpr_log = np.log10(np.clip(fpr, lam, 1.0))
    
    # Normalized integral between log10(lam) and log10(1) = 0
    total_log_span = -np.log10(lam)
    log_auc = _trapezoid_integral(tpr, fpr_log) / total_log_span
    return float(np.clip(log_auc, 0.0, 1.0))


def calculate_optimal_mcc(
    labels: np.ndarray,
    scores: np.ndarray,
    lower_is_better: bool = True
) -> Tuple[float, float]:
    """
    Computes the maximum Matthews Correlation Coefficient (MCC) across all candidate score thresholds.

    Parameters
    ----------
    labels : np.ndarray
        Binary array.
    scores : np.ndarray
        Docking scores.
    lower_is_better : bool, default True
        Sorting order.

    Returns
    -------
    best_mcc : float
        Optimal MCC in [-1, +1].
    best_threshold : float
        Docking score threshold achieving optimal MCC.
    """
    y_true = np.asarray(labels, dtype=int)
    y_scores = -np.asarray(scores, dtype=float) if lower_is_better else np.asarray(scores, dtype=float)
    
    if len(np.unique(y_true)) < 2:
        return 0.0, 0.0
        
    fpr, tpr, thresholds = metrics.roc_curve(y_true, y_scores)
    
    n_pos = np.sum(y_true == 1)
    n_neg = np.sum(y_true == 0)
    
    best_mcc = -1.0
    best_thresh = 0.0
    
    for i in range(len(thresholds)):
        tp = tpr[i] * n_pos
        fp = fpr[i] * n_neg
        fn = n_pos - tp
        tn = n_neg - fp
        
        denom = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
        if denom == 0:
            mcc = 0.0
        else:
            mcc = (tp * tn - fp * fn) / denom
            
        if mcc > best_mcc:
            best_mcc = mcc
            # Revert negation if lower_is_better
            best_thresh = -thresholds[i] if lower_is_better else thresholds[i]
            
    return float(best_mcc), float(best_thresh)


def evaluate_all_enrichment_metrics(
    labels: np.ndarray,
    scores: np.ndarray,
    lower_is_better: bool = True
) -> Dict[str, Any]:
    """
    Evaluates comprehensive virtual screening enrichment statistics.

    Returns
    -------
    metrics_dict : dict
        Consolidated dictionary of all enrichment metrics.
    """
    roc_auc = calculate_roc_auc(labels, scores, lower_is_better)
    pr_auc = calculate_pr_auc(labels, scores, lower_is_better)
    bedroc_20 = calculate_bedroc(labels, scores, alpha=20.0, lower_is_better=lower_is_better)
    bedroc_80 = calculate_bedroc(labels, scores, alpha=80.5, lower_is_better=lower_is_better)
    rie_20 = calculate_rie(labels, scores, alpha=20.0, lower_is_better=lower_is_better)
    
    ef_1, ef_1_max = calculate_enrichment_factor(labels, scores, fraction=0.01, lower_is_better=lower_is_better)
    ef_5, ef_5_max = calculate_enrichment_factor(labels, scores, fraction=0.05, lower_is_better=lower_is_better)
    ef_10, ef_10_max = calculate_enrichment_factor(labels, scores, fraction=0.10, lower_is_better=lower_is_better)
    
    log_auc = calculate_log_auc(labels, scores, lower_is_better=lower_is_better)
    mcc, mcc_thresh = calculate_optimal_mcc(labels, scores, lower_is_better=lower_is_better)
    
    n_actives = int(np.sum(labels))
    n_total = len(labels)
    n_decoys = n_total - n_actives
    
    return {
        "n_total": n_total,
        "n_actives": n_actives,
        "n_decoys": n_decoys,
        "roc_auc": roc_auc,
        "pr_auc": pr_auc,
        "bedroc_20": bedroc_20,
        "bedroc_80": bedroc_80,
        "rie_20": rie_20,
        "ef_1pct": ef_1,
        "ef_1pct_max": ef_1_max,
        "ef_5pct": ef_5,
        "ef_5pct_max": ef_5_max,
        "ef_10pct": ef_10,
        "ef_10pct_max": ef_10_max,
        "log_auc": log_auc,
        "optimal_mcc": mcc,
        "optimal_mcc_threshold": mcc_thresh
    }
