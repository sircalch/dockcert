"""
Heavy-atom and symmetry-corrected Root-Mean-Square Deviation (RMSD) for ligand poses.
"""

from typing import Dict, Any, Optional, List
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
        Docked pose coordinates (shape: [N, 3]), same atom order as the reference.

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


_COV_RADII = {"C": 0.76, "N": 0.71, "O": 0.66, "S": 1.05, "P": 1.07, "F": 0.57, "CL": 1.02,
              "BR": 1.20, "I": 1.39, "B": 0.84, "SE": 1.20, "SI": 1.11}


def infer_bonds(coords: np.ndarray, elements: List[str], tolerance: float = 0.45) -> List[set]:
    """Adjacency lists from interatomic distances (covalent radii + tolerance, in Angstrom)."""
    c = np.asarray(coords, dtype=np.float64)
    r = np.array([_COV_RADII.get(e.upper(), 0.77) for e in elements])
    d = np.linalg.norm(c[:, None, :] - c[None, :, :], axis=-1)
    bonded = (d < r[:, None] + r[None, :] + tolerance) & (d > 0.4)
    np.fill_diagonal(bonded, False)
    return [set(np.nonzero(row)[0].tolist()) for row in bonded]


def _isomorphisms(el_a, adj_a, el_b, adj_b, limit):
    """Yields element- and bond-preserving maps a -> b (backtracking over a in BFS order)."""
    n = len(el_a)
    order, seen = [], set()
    for root in range(n):
        if root in seen:
            continue
        queue = [root]
        seen.add(root)
        while queue:
            i = queue.pop(0)
            order.append(i)
            for j in sorted(adj_a[i]):
                if j not in seen:
                    seen.add(j)
                    queue.append(j)
    mapping, used = {}, set()
    count = [0]

    def extend(k):
        if count[0] >= limit:
            return
        if k == n:
            count[0] += 1
            yield dict(mapping)
            return
        i = order[k]
        mapped_nb = [j for j in adj_a[i] if j in mapping]
        for cand in range(n):
            if cand in used or el_b[cand] != el_a[i] or len(adj_b[cand]) != len(adj_a[i]):
                continue
            if any(mapping[j] not in adj_b[cand] for j in mapped_nb):
                continue
            mapping[i] = cand
            used.add(cand)
            yield from extend(k + 1)
            used.discard(cand)
            del mapping[i]

    yield from extend(0)


def calculate_symmetry_corrected_rmsd(
    coords_ref: np.ndarray,
    coords_dock: np.ndarray,
    elements: Optional[List[str]] = None,
    elements_dock: Optional[List[str]] = None,
    max_mappings: int = 200000,
) -> float:
    """
    Symmetry-corrected heavy-atom RMSD without superposition (the redocking convention).

    Both poses are turned into molecular graphs (bonds inferred from distances) and the RMSD is
    minimised over all element- and bond-preserving atom correspondences (graph isomorphisms).
    This handles symmetric groups and a different atom order between the crystal ligand and the
    docking output, and matches the definition of RDKit CalcRMS and spyrmsd.

    If the two graphs are not isomorphic (e.g. a missing atom or a distorted bond), it falls back
    to a per-element Hungarian assignment, which ignores bonding and is only a LOWER BOUND.

    Parameters
    ----------
    coords_ref, coords_dock : np.ndarray
        Heavy-atom coordinates (shape [N, 3]); the atom order may differ.
    elements : list of str, optional
        Element symbols of the reference atoms. Without it, the plain RMSD is returned.
    elements_dock : list of str, optional
        Element symbols of the docked atoms (default: the reference list, i.e. same atom order).
    max_mappings : int
        Cap on the number of enumerated isomorphisms.

    Returns
    -------
    rmsd_sym : float
        Symmetry-corrected RMSD in Angstroms.
    """
    c_ref = np.asarray(coords_ref, dtype=np.float64)
    c_dock = np.asarray(coords_dock, dtype=np.float64)
    if c_ref.shape != c_dock.shape or len(c_ref) == 0:
        raise ValueError(f"Heavy-atom count mismatch between reference and pose: {c_ref.shape} vs {c_dock.shape}")
    n_atoms = len(c_ref)
    if elements is None or len(elements) != n_atoms:
        return calculate_heavy_atom_rmsd(c_ref, c_dock)
    el_ref = [e.upper() for e in elements]
    el_dock = [e.upper() for e in (elements_dock if elements_dock is not None else elements)]
    if sorted(el_ref) != sorted(el_dock):
        raise ValueError("Reference and pose have different heavy-atom compositions")

    adj_ref = infer_bonds(c_ref, el_ref)
    adj_dock = infer_bonds(c_dock, el_dock)
    best = None
    for m in _isomorphisms(el_ref, adj_ref, el_dock, adj_dock, max_mappings):
        idx = [m[i] for i in range(n_atoms)]
        sq = float(np.sum((c_ref - c_dock[idx]) ** 2))
        if best is None or sq < best:
            best = sq
    if best is not None:
        return float(np.sqrt(best / n_atoms))

    # Fallback: per-element Hungarian assignment (lower bound; bonding not enforced)
    total_sq_dist = 0.0
    for elem in set(el_ref):
        ia = [i for i, e in enumerate(el_ref) if e == elem]
        ib = [i for i, e in enumerate(el_dock) if e == elem]
        cost = np.sum((c_ref[ia][:, None, :] - c_dock[ib][None, :, :]) ** 2, axis=-1)
        r, c = linear_sum_assignment(cost)
        total_sq_dist += float(np.sum(cost[r, c]))
    return float(np.sqrt(total_sq_dist / n_atoms))


def evaluate_redocking_success(
    rmsd_values: List[float],
    threshold_pass: float = 2.0,
    threshold_warn: float = 3.0,
    ranked_poses: bool = False,
) -> Dict[str, Any]:
    """
    Evaluates redocking against the 2 A pose-recovery criterion.

    * ranked_poses=True: rmsd_values are the poses of ONE docking run in score order. The verdict
      uses the top-scored pose, which is the pose a user of the result would take. The best of the
      N poses is reported separately because it measures sampling, not scoring.
      PASS: top-1 <= 2 A. WARNING: a lower-ranked pose <= 2 A (scoring failure) or top-1 <= 3 A.
      FAIL otherwise.
    * ranked_poses=False: rmsd_values are the top-scored poses of a benchmark of complexes.
      PASS: success rate (<= 2 A) >= 70 %. WARNING: >= 50 %. FAIL otherwise.

    Returns
    -------
    result : dict
        n_poses, top_rmsd, mean/median/min RMSD, success_rate_2a, status and recommendation.
    """
    rmsds = np.asarray(rmsd_values, dtype=np.float64)
    n = len(rmsds)
    if n == 0:
        return {
            "n_poses": 0,
            "top_rmsd": float("nan"),
            "mean_rmsd": float("nan"),
            "median_rmsd": float("nan"),
            "min_rmsd": float("nan"),
            "success_rate_2a": 0.0,
            "status": "FAIL",
            "recommendation": "No poses provided for redocking validation."
        }

    top = float(rmsds[0])
    min_rmsd = float(np.min(rmsds))
    n_pass = int(np.sum(rmsds <= threshold_pass))
    success_rate = n_pass / n * 100.0

    if ranked_poses:
        if top <= threshold_pass:
            status = "PASS"
            recommendation = f"Top-scored pose reproduces the crystallographic binding mode (RMSD = {top:.2f} A <= {threshold_pass} A)."
        elif min_rmsd <= threshold_pass:
            status = "WARNING"
            recommendation = (f"Sampling found the native mode (best of {n} poses: {min_rmsd:.2f} A) but the scoring "
                              f"function ranked a wrong pose first (top-1 RMSD = {top:.2f} A).")
        elif top <= threshold_warn:
            status = "WARNING"
            recommendation = f"Borderline redocking accuracy (top-1 RMSD = {top:.2f} A). Inspect the box, protonation and ligand preparation."
        else:
            status = "FAIL"
            recommendation = f"Redocking failed to reproduce the crystallographic binding mode (top-1 RMSD = {top:.2f} A, best of {n} = {min_rmsd:.2f} A)."
    else:
        if success_rate >= 70.0:
            status = "PASS"
        elif success_rate >= 50.0:
            status = "WARNING"
        else:
            status = "FAIL"
        recommendation = f"Redocking success rate {success_rate:.1f}% ({n_pass}/{n} complexes with top-1 RMSD <= {threshold_pass} A)."

    return {
        "n_poses": n,
        "top_rmsd": top,
        "mean_rmsd": float(np.mean(rmsds)),
        "median_rmsd": float(np.median(rmsds)),
        "min_rmsd": min_rmsd,
        "success_rate_2a": success_rate,
        "status": status,
        "recommendation": recommendation
    }
