"""
Decoy bias and physicochemical property distribution diagnostics.
"""

from typing import Dict, Any, Optional
import numpy as np
from scipy import stats


def evaluate_decoy_bias(
    active_properties: Dict[str, np.ndarray],
    decoy_properties: Dict[str, np.ndarray]
) -> Dict[str, Any]:
    """
    Evaluates potential artificial enrichment bias by comparing property distributions
    (e.g., MW, LogP, HBD, HBA, Rotatable Bonds) between actives and decoys.

    Parameters
    ----------
    active_properties : dict
        Mapping of property name -> 1D array for active molecules.
    decoy_properties : dict
        Mapping of property name -> 1D array for decoy molecules.

    Returns
    -------
    result : dict
        Property-by-property Kolmogorov-Smirnov and Wasserstein statistics,
        bias risk classification (LOW / MODERATE / HIGH), and recommendations.
    """
    common_props = set(active_properties.keys()).intersection(set(decoy_properties.keys()))
    if not common_props:
        return {
            "status": "PASS",
            "risk_level": "LOW",
            "properties_evaluated": 0,
            "property_metrics": {},
            "recommendation": "No property metadata provided for decoy bias auditing. Ensure decoys are property-matched."
        }
        
    property_metrics = {}
    ks_pvalues = []
    
    for prop in common_props:
        a_vals = np.asarray(active_properties[prop], dtype=float)
        d_vals = np.asarray(decoy_properties[prop], dtype=float)
        
        if len(a_vals) < 2 or len(d_vals) < 2:
            continue
            
        ks_res = stats.ks_2samp(a_vals, d_vals)
        w_dist = stats.wasserstein_distance(a_vals, d_vals)
        
        property_metrics[prop] = {
            "active_mean": float(np.mean(a_vals)),
            "decoy_mean": float(np.mean(d_vals)),
            "ks_statistic": float(ks_res.statistic),
            "ks_pvalue": float(ks_res.pvalue),
            "wasserstein_distance": float(w_dist)
        }
        ks_pvalues.append(ks_res.pvalue)

    # If KS p-value is extremely low (< 1e-4) across multiple physical properties, decoys are poorly matched
    significant_biases = sum(1 for p in ks_pvalues if p < 0.001)
    
    if significant_biases == 0:
        risk_level = "LOW"
        status = "PASS"
        recommendation = "Decoys and actives are well property-matched. Low risk of artificial enrichment bias."
    elif significant_biases <= 2:
        risk_level = "MODERATE"
        status = "WARNING"
        recommendation = "Moderate distributional discrepancy observed in physicochemical properties between actives and decoys."
    else:
        risk_level = "HIGH"
        status = "WARNING"
        recommendation = "High risk of decoy bias: actives and decoys differ significantly in physical properties (e.g. MW, LogP). Enrichment may be artificially inflated."

    return {
        "status": status,
        "risk_level": risk_level,
        "properties_evaluated": len(property_metrics),
        "property_metrics": property_metrics,
        "recommendation": recommendation
    }
