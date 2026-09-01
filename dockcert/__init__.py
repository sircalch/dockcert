"""
DockCert: Automated Statistical Validation, Enrichment Metrics, and Reproducibility
Assessment for Molecular Docking and Virtual Screening Studies.
"""

__version__ = "1.0.0"
__author__ = "Andres Monreal-Hernández"
__license__ = "MIT"

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
from dockcert.core.rmsd import (
    calculate_heavy_atom_rmsd,
    calculate_symmetry_corrected_rmsd,
    evaluate_redocking_success
)
from dockcert.core.bias import evaluate_decoy_bias
from dockcert.core.bootstrap import bootstrap_enrichment_ci
from dockcert.core.scoring import assess_docking_quality, DockingValidationReport

__all__ = [
    "__version__",
    "calculate_roc_auc",
    "calculate_pr_auc",
    "calculate_bedroc",
    "calculate_rie",
    "calculate_enrichment_factor",
    "calculate_log_auc",
    "calculate_optimal_mcc",
    "evaluate_all_enrichment_metrics",
    "calculate_heavy_atom_rmsd",
    "calculate_symmetry_corrected_rmsd",
    "evaluate_redocking_success",
    "evaluate_decoy_bias",
    "bootstrap_enrichment_ci",
    "assess_docking_quality",
    "DockingValidationReport"
]
