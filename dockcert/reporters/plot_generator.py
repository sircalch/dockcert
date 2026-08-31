"""
Publication-ready vector figure generation for docking and virtual screening validation.
"""

from typing import List, Optional
import os
import numpy as np
import matplotlib.pyplot as plt
from sklearn import metrics
from dockcert.core.scoring import DockingValidationReport

# Clean scientific aesthetics
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.size': 11,
    'axes.labelsize': 12,
    'axes.titlesize': 13,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 10,
    'figure.titlesize': 14,
    'figure.dpi': 300,
    'lines.linewidth': 2.0,
    'grid.alpha': 0.3,
    'grid.linestyle': '--'
})


def generate_docking_figures(
    labels: np.ndarray,
    scores: np.ndarray,
    report: DockingValidationReport,
    output_dir: str,
    lower_is_better: bool = True,
    rmsd_values: Optional[List[float]] = None,
    formats: List[str] = ("png", "svg", "pdf")
) -> List[str]:
    """
    Generates high-resolution publication figures for virtual screening and docking validation.

    Parameters
    ----------
    labels : np.ndarray
        Binary labels (1 = Active, 0 = Decoy).
    scores : np.ndarray
        Docking scores.
    report : DockingValidationReport
        Validation report.
    output_dir : str
        Target output directory.
    lower_is_better : bool
        Sort order.
    rmsd_values : list of float, optional
        Pose RMSD values.
    formats : list of str
        Formats to save.

    Returns
    -------
    saved_paths : list of str
        List of generated file paths.
    """
    os.makedirs(output_dir, exist_ok=True)
    saved_files = []
    
    y_true = np.asarray(labels, dtype=int)
    y_scores = -np.asarray(scores, dtype=float) if lower_is_better else np.asarray(scores, dtype=float)
    
    # 1. Main 4-Panel Validation Panel
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # Subplot (0, 0): ROC Curve
    ax_roc = axes[0, 0]
    fpr, tpr, _ = metrics.roc_curve(y_true, y_scores)
    roc_item = report.enrichment_metrics.get("roc_auc")
    roc_val = roc_item.value if roc_item else metrics.roc_auc_score(y_true, y_scores)
    ci_str = f" [{roc_item.ci_lower_95:.2f}, {roc_item.ci_upper_95:.2f}]" if roc_item and roc_item.ci_lower_95 else ""
    
    ax_roc.plot(fpr, tpr, color="#2563eb", label=f"Docking Model (AUC = {roc_val:.3f}{ci_str})")
    ax_roc.plot([0, 1], [0, 1], color="#94a3b8", linestyle="--", label="Random Baseline (AUC = 0.500)")
    ax_roc.set_xlabel("False Positive Rate (1 - Specificity)")
    ax_roc.set_ylabel("True Positive Rate (Sensitivity)")
    ax_roc.set_title("Receiver Operating Characteristic (ROC)")
    ax_roc.grid(True)
    ax_roc.legend(loc="lower right", frameon=True)
    
    # Subplot (0, 1): Precision-Recall Curve
    ax_pr = axes[0, 1]
    precision, recall, _ = metrics.precision_recall_curve(y_true, y_scores)
    pr_item = report.enrichment_metrics.get("pr_auc")
    pr_val = pr_item.value if pr_item else metrics.average_precision_score(y_true, y_scores)
    baseline_pr = float(np.mean(y_true))
    
    ax_pr.plot(recall, precision, color="#10b981", label=f"Precision-Recall (PR-AUC = {pr_val:.3f})")
    ax_pr.axhline(baseline_pr, color="#94a3b8", linestyle="--", label=f"Active Proportion ({baseline_pr*100:.1f}%)")
    ax_pr.set_xlabel("Recall (True Positive Rate)")
    ax_pr.set_ylabel("Precision (Positive Predictive Value)")
    ax_pr.set_title("Precision-Recall (PR) Curve")
    ax_pr.grid(True)
    ax_pr.legend(loc="upper right", frameon=True)

    # Subplot (1, 0): Score Distribution Separation
    ax_dist = axes[1, 0]
    active_scores = scores[y_true == 1]
    decoy_scores = scores[y_true == 0]
    
    ax_dist.hist(decoy_scores, bins=30, density=True, alpha=0.5, color="#64748b", label=f"Decoys (N = {len(decoy_scores)})")
    ax_dist.hist(active_scores, bins=30, density=True, alpha=0.6, color="#e11d48", label=f"Actives (N = {len(active_scores)})")
    ax_dist.set_xlabel("Docking Score / Binding Affinity (kcal/mol)")
    ax_dist.set_ylabel("Probability Density")
    ax_dist.set_title("Score Distribution Separation")
    ax_dist.grid(True)
    ax_dist.legend(loc="upper right", frameon=True)

    # Subplot (1, 1): Active Recovery Rate vs Database Screened
    ax_rec = axes[1, 1]
    order = np.argsort(scores) if lower_is_better else np.argsort(-scores)
    sorted_labels = y_true[order]
    cum_actives = np.cumsum(sorted_labels) / np.sum(sorted_labels) * 100.0
    db_fraction = (np.arange(len(sorted_labels)) + 1) / len(sorted_labels) * 100.0
    
    ax_rec.plot(db_fraction, cum_actives, color="#8b5cf6", label="Observed Active Recovery")
    ax_rec.plot([0, 100], [0, 100], color="#94a3b8", linestyle="--", label="Random Selection")
    ax_rec.axvline(1.0, color="#f59e0b", linestyle=":", label="Top 1% Cutoff")
    ax_rec.axvline(5.0, color="#ec4899", linestyle=":", label="Top 5% Cutoff")
    ax_rec.set_xlabel("Database Screened (%)")
    ax_rec.set_ylabel("Actives Recovered (%)")
    ax_rec.set_title("Virtual Screening Cumulative Enrichment")
    ax_rec.grid(True)
    ax_rec.legend(loc="lower right", frameon=True)

    plt.tight_layout()
    for fmt in formats:
        p = os.path.join(output_dir, f"dockcert_validation_overview.{fmt}")
        plt.savefig(p, dpi=300, bbox_inches="tight")
        saved_files.append(p)
    plt.close()

    # 2. Redocking RMSD Figure (if provided)
    if rmsd_values is not None and len(rmsd_values) > 0:
        fig, ax = plt.subplots(figsize=(7, 4.5))
        rmsds = np.asarray(rmsd_values)
        ax.hist(rmsds, bins=20, color="#0284c7", edgecolor="black", alpha=0.7)
        ax.axvline(2.0, color="#dc2626", linestyle="--", linewidth=2.0, label="2.0 A Success Threshold")
        ax.set_xlabel("Pose RMSD (Angstroms)")
        ax.set_ylabel("Frequency")
        ax.set_title(f"Redocking Pose Accuracy (Best RMSD = {np.min(rmsds):.2f} A)")
        ax.grid(True)
        ax.legend(frameon=True)
        plt.tight_layout()
        
        for fmt in formats:
            p = os.path.join(output_dir, f"dockcert_redocking_rmsd.{fmt}")
            plt.savefig(p, dpi=300, bbox_inches="tight")
            saved_files.append(p)
        plt.close()

    return saved_files
