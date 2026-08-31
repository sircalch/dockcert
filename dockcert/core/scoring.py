"""
Docking and Virtual Screening Validation Scoring Engine.
"""

from typing import Dict, Any, Optional, List
from dataclasses import dataclass, asdict
import numpy as np

from dockcert.core.enrichment import evaluate_all_enrichment_metrics
from dockcert.core.rmsd import evaluate_redocking_success
from dockcert.core.bias import evaluate_decoy_bias
from dockcert.core.bootstrap import bootstrap_enrichment_ci


@dataclass
class MetricValidationItem:
    name: str
    value: float
    ci_lower_95: Optional[float]
    ci_upper_95: Optional[float]
    threshold_pass: float
    threshold_warn: float
    status: str
    message: str


@dataclass
class DockingValidationReport:
    overall_status: str
    validation_score: str
    redocking_result: Optional[Dict[str, Any]]
    enrichment_metrics: Dict[str, MetricValidationItem]
    decoy_bias_result: Optional[Dict[str, Any]]
    dataset_summary: Dict[str, Any]
    recommendations: List[str]
    provenance: Dict[str, Any]

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


def assess_docking_quality(
    labels: Optional[np.ndarray] = None,
    scores: Optional[np.ndarray] = None,
    rmsd_values: Optional[List[float]] = None,
    active_properties: Optional[Dict[str, np.ndarray]] = None,
    decoy_properties: Optional[Dict[str, np.ndarray]] = None,
    lower_is_better: bool = True,
    roc_auc_pass: float = 0.70,
    roc_auc_warn: float = 0.60,
    bedroc_pass: float = 0.50,
    bedroc_warn: float = 0.30,
    ef1_pass: float = 5.0,
    ef1_warn: float = 2.0
) -> DockingValidationReport:
    """
    Evaluates comprehensive quality assurance for molecular docking and virtual screening.

    Parameters
    ----------
    labels : np.ndarray, optional
        Binary array of active/decoy classifications.
    scores : np.ndarray, optional
        Docking scores.
    rmsd_values : list of float, optional
        Pose RMSD values from redocking/cross-docking experiments.
    active_properties : dict, optional
        Properties for active ligands.
    decoy_properties : dict, optional
        Properties for decoy ligands.
    lower_is_better : bool, default True
        Sort order.

    Returns
    -------
    report : DockingValidationReport
        Consolidated validation report.
    """
    statuses = []
    recommendations = []
    
    # 1. Redocking Evaluation
    redocking_res = None
    if rmsd_values is not None and len(rmsd_values) > 0:
        redocking_res = evaluate_redocking_success(rmsd_values)
        statuses.append(redocking_res["status"])
        if redocking_res["status"] != "PASS":
            recommendations.append(redocking_res["recommendation"])

    # 2. Enrichment Evaluation
    enrichment_items: Dict[str, MetricValidationItem] = {}
    dataset_summary: Dict[str, Any] = {}
    
    if labels is not None and scores is not None and len(labels) > 0:
        y_true = np.asarray(labels, dtype=int)
        y_scores = np.asarray(scores, dtype=float)
        
        raw_metrics = evaluate_all_enrichment_metrics(y_true, y_scores, lower_is_better=lower_is_better)
        ci_results = bootstrap_enrichment_ci(y_true, y_scores, lower_is_better=lower_is_better)
        
        dataset_summary = {
            "n_total": raw_metrics["n_total"],
            "n_actives": raw_metrics["n_actives"],
            "n_decoys": raw_metrics["n_decoys"],
            "active_ratio_pct": float(raw_metrics["n_actives"] / raw_metrics["n_total"] * 100.0)
        }
        
        # Helper to register and test metrics
        def check_metric(key, name, val, th_pass, th_warn, greater_is_better=True):
            ci_low, ci_high = None, None
            if key in ci_results:
                _, ci_low, ci_high = ci_results[key]
                
            if greater_is_better:
                if val >= th_pass:
                    st = "PASS"
                    msg = f"Strong discrimination ({name} = {val:.2f} >= {th_pass:.2f})."
                elif val >= th_warn:
                    st = "WARNING"
                    msg = f"Moderate discrimination ({name} = {val:.2f} >= {th_warn:.2f})."
                else:
                    st = "FAIL"
                    msg = f"Inadequate discrimination ({name} = {val:.2f} < {th_warn:.2f})."
            else:
                if val <= th_pass:
                    st = "PASS"
                    msg = f"Optimal value ({name} = {val:.2f} <= {th_pass:.2f})."
                elif val <= th_warn:
                    st = "WARNING"
                    msg = f"Borderline value ({name} = {val:.2f})."
                else:
                    st = "FAIL"
                    msg = f"Exceeds acceptable threshold ({name} = {val:.2f} > {th_warn:.2f})."
                    
            statuses.append(st)
            enrichment_items[key] = MetricValidationItem(
                name=name,
                value=float(val),
                ci_lower_95=ci_low,
                ci_upper_95=ci_high,
                threshold_pass=float(th_pass),
                threshold_warn=float(th_warn),
                status=st,
                message=msg
            )

        check_metric("roc_auc", "ROC-AUC", raw_metrics["roc_auc"], roc_auc_pass, roc_auc_warn)
        check_metric("bedroc_20", "BEDROC (alpha=20.0)", raw_metrics["bedroc_20"], bedroc_pass, bedroc_warn)
        check_metric("ef_1pct", "EF 1%", raw_metrics["ef_1pct"], ef1_pass, ef1_warn)
        check_metric("ef_5pct", "EF 5%", raw_metrics["ef_5pct"], 3.0, 1.5)
        check_metric("optimal_mcc", "Optimal MCC", raw_metrics["optimal_mcc"], 0.40, 0.20)

    # 3. Decoy Bias Audit
    bias_res = None
    if active_properties is not None and decoy_properties is not None:
        bias_res = evaluate_decoy_bias(active_properties, decoy_properties)
        statuses.append(bias_res["status"])
        if bias_res["status"] != "PASS":
            recommendations.append(bias_res["recommendation"])

    # Overall scoring decision
    if "FAIL" in statuses:
        overall_status = "FAIL"
        validation_score = "VALIDATION STATUS = UNVERIFIED / REJECTED"
    elif "WARNING" in statuses:
        overall_status = "WARNING"
        validation_score = "VALIDATION STATUS = BORDERLINE / ACCEPTABLE WITH LIMITATIONS"
    else:
        overall_status = "PASS"
        validation_score = "VALIDATION STATUS = FULLY CERTIFIED"

    return DockingValidationReport(
        overall_status=overall_status,
        validation_score=validation_score,
        redocking_result=redocking_res,
        enrichment_metrics=enrichment_items,
        decoy_bias_result=bias_res,
        dataset_summary=dataset_summary,
        recommendations=recommendations,
        provenance={
            "tool": "DockCert",
            "version": "1.0.0",
            "citation": "Monreal-Hernández, A. (2026). DockCert: An Open-Source Toolkit for Statistical Validation, Enrichment Metrics, and Reproducibility Assessment of Molecular Docking Studies."
        }
    )
