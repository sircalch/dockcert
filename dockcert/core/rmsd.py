"""
Heavy-atom and symmetry-corrected Root-Mean-Square Deviation (RMSD) for ligand poses.
"""

from typing import Dict, Any, Optional, List, Tuple
import numpy as np
from scipy.optimize import linear_sum_assignment


def calculate_heavy_atom_rmsd(
    coords_ref: np.ndarray,
    coords_dock: np.ndarray
) -> float:
    """
    Computes standard Cartesian Root-Mean-Square Deviation (RMSD) between two matching atom sets.

    Parameters
    ----------
    coords_ref : np.ndarray
        Reference coordinates (shape: [N, 3]).
    coords_dock : np.ndarray
        Docked pose coordinates (shape: [N, 3]).

    Returns
    -------
    rmsd : float
        RMSD in Angstroms.
    """
    c_ref = np.asarray(coords_ref, dtype=np.float64)
    c_dock = np.asarray(coords_dock, dtype=np.float64)
    
    if c_ref.shape != c_dock.shape or len(c_ref) == 0:
        raise ValueError(f"Coordinate shape mismatch: {c_ref.shape} vs {c_dock.shape}")
        
    diff = c_ref - c_dock
    sq_dist = np.sum(diff ** 2, axis=-1)
    mean_sq = np.mean(sq_dist)
    return float(np.sqrt(mean_sq))


def calculate_symmetry_corrected_rmsd(
    coords_ref: np.ndarray,
    coords_dock: np.ndarray,
    elements: Optional[List[str]] = None
) -> float:
    """
    Computes symmetry-corrected RMSD by solving the minimum-weight bipartite matching
    (Hungarian algorithm) across symmetric / identical chemical elements.

    Parameters
    ----------
    coords_ref : np.ndarray
        Reference coordinates (shape: [N, 3]).
    coords_dock : np.ndarray
        Docked pose coordinates (shape: [N, 3]).
    elements : list of str, optional
        Element symbols for each atom (e.g. ['C', 'C', 'O', 'N']).

    Returns
    -------
    rmsd_sym : float
        Symmetry-corrected RMSD in Angstroms.
    """
    c_ref = np.asarray(coords_ref, dtype=np.float64)
    c_dock = np.asarray(coords_dock, dtype=np.float64)
    n_atoms = len(c_ref)
    
    if elements is None or len(elements) != n_atoms:
        # Default standard RMSD
        return calculate_heavy_atom_rmsd(c_ref, c_dock)
        
    unique_elems = np.unique(elements)
    total_sq_dist = 0.0
    
    for elem in unique_elems:
        indices = [i for i, el in enumerate(elements) if el == elem]
        if not indices:
            continue
            
        sub_ref = c_ref[indices]
        sub_dock = c_dock[indices]
        
        # Distance cost matrix between ref atoms and dock atoms of this element
        # cost[i, j] = ||sub_ref[i] - sub_dock[j]||^2
        cost = np.sum((sub_ref[:, np.newaxis, :] - sub_dock[np.newaxis, :, :]) ** 2, axis=-1)
        
        row_ind, col_ind = linear_sum_assignment(cost)
        total_sq_dist += np.sum(cost[row_ind, col_ind])
        
    rmsd_sym = np.sqrt(total_sq_dist / n_atoms)
    return float(rmsd_sym)


def evaluate_redocking_success(
    rmsd_values: List[float],
    threshold_pass: float = 2.0,
    threshold_warn: float = 3.0
) -> Dict[str, Any]:
    """
    Evaluates redocking success rate against established CADD validation standards.

    Parameters
    ----------
    rmsd_values : list of float
        RMSD values across poses or benchmark targets.
    threshold_pass : float, default 2.0
        Standard threshold for successful pose recovery (2.0 A).
    threshold_warn : float, default 3.0
        Threshold for borderline pose prediction (3.0 A).

    Returns
    -------
    result : dict
        Mean RMSD, median RMSD, success rate percentage, and PASS/WARNING/FAIL status.
    """
    rmsds = np.asarray(rmsd_values, dtype=np.float64)
    n = len(rmsds)
    if n == 0:
        return {
            "n_poses": 0,
            "mean_rmsd": 0.0,
            "median_rmsd": 0.0,
            "min_rmsd": 0.0,
            "success_rate_2a": 0.0,
            "status": "FAIL",
            "recommendation": "No poses provided for redocking validation."
        }
        
    mean_rmsd = float(np.mean(rmsds))
    median_rmsd = float(np.median(rmsds))
    min_rmsd = float(np.min(rmsds))
    
    n_pass = int(np.sum(rmsds <= threshold_pass))
    success_rate = (n_pass / n) * 100.0
    
    # Evaluate status
    if min_rmsd <= threshold_pass:
        if success_rate >= 80.0 or min_rmsd <= 1.5:
            status = "PASS"
            recommendation = f"Excellent redocking accuracy (Best RMSD = {min_rmsd:.2f} A <= 2.0 A)."
        else:
            status = "WARNING"
            recommendation = f"Best pose reproduced native binding mode (RMSD = {min_rmsd:.2f} A), but ensemble success rate is moderate ({success_rate:.1f}%)."
    elif min_rmsd <= threshold_warn:
        status = "WARNING"
        recommendation = f"Borderline redocking accuracy (Best RMSD = {min_rmsd:.2f} A). Consider inspecting scoring function or box definition."
    else:
        status = "FAIL"
        recommendation = f"Redocking failed to reproduce experimental binding mode (Best RMSD = {min_rmsd:.2f} A > 3.0 A)."

    return {
        "n_poses": n,
        "mean_rmsd": mean_rmsd,
        "median_rmsd": median_rmsd,
        "min_rmsd": min_rmsd,
        "success_rate_2a": success_rate,
        "status": status,
        "recommendation": recommendation
    }
