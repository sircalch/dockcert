"""
Tests for Cartesian and symmetry-corrected RMSD calculations.
"""

import numpy as np
import pytest
from dockcert.core.rmsd import (
    calculate_heavy_atom_rmsd,
    calculate_symmetry_corrected_rmsd,
    evaluate_redocking_success
)


def test_heavy_atom_rmsd_identity():
    coords = np.array([
        [0.0, 0.0, 0.0],
        [1.5, 0.0, 0.0],
        [1.5, 1.5, 0.0]
    ])
    rmsd = calculate_heavy_atom_rmsd(coords, coords)
    assert np.isclose(rmsd, 0.0)
    
    # 1.0 A translation along Z
    shifted = coords + np.array([0.0, 0.0, 1.0])
    rmsd_shift = calculate_heavy_atom_rmsd(coords, shifted)
    assert np.isclose(rmsd_shift, 1.0)


def test_symmetry_corrected_rmsd():
    # Symmetric 2-oxygen group (e.g. carboxylate)
    # Ref: O1 at (0, 1, 0), O2 at (0, -1, 0)
    # Dock: O1 at (0, -1, 0), O2 at (0, 1, 0) (flipped indices)
    c_ref = np.array([[0.0, 1.0, 0.0], [0.0, -1.0, 0.0]])
    c_dock = np.array([[0.0, -1.0, 0.0], [0.0, 1.0, 0.0]])
    elements = ["O", "O"]
    
    # Standard RMSD without symmetry would be 2.0 A
    std_rmsd = calculate_heavy_atom_rmsd(c_ref, c_dock)
    assert np.isclose(std_rmsd, 2.0)
    
    # Symmetry-corrected RMSD matches equivalent oxygens -> 0.0 A
    sym_rmsd = calculate_symmetry_corrected_rmsd(c_ref, c_dock, elements=elements)
    assert np.isclose(sym_rmsd, 0.0)


def test_evaluate_redocking_success_benchmark():
    """Top-1 RMSDs across complexes: verdict on the success rate, not on the single best value."""
    res = evaluate_redocking_success([1.2, 1.4, 1.1, 0.9, 1.8, 2.5, 3.2])
    assert res["status"] == "PASS"            # 5/7 = 71 %
    assert res["min_rmsd"] == 0.9
    assert evaluate_redocking_success([1.2, 1.4, 2.5, 3.2])["status"] == "WARNING"   # 50 %
    assert evaluate_redocking_success([1.2, 4.5, 5.2, 6.1])["status"] == "FAIL"     # one lucky target is not a pass


def test_evaluate_redocking_ranked_poses():
    assert evaluate_redocking_success([0.4, 3.4, 11.0], ranked_poses=True)["status"] == "PASS"
    res = evaluate_redocking_success([4.4, 1.5, 10.0], ranked_poses=True)
    assert res["status"] == "WARNING" and res["top_rmsd"] == 4.4 and res["min_rmsd"] == 1.5
    assert evaluate_redocking_success([4.4, 10.2, 10.8], ranked_poses=True)["status"] == "FAIL"


def test_graph_rmsd_is_order_independent():
    """Benzene-like ring with shuffled atom order and a 180-degree flip: RMSD must be 0."""
    ang = np.linspace(0, 2 * np.pi, 7)[:-1]
    ring = np.c_[1.39 * np.cos(ang), 1.39 * np.sin(ang), np.zeros(6)]
    ring = np.vstack([ring, [[2.9, 0.0, 0.0]]])          # substituent on atom 0
    el = ["C"] * 6 + ["N"]
    perm = [3, 1, 5, 0, 4, 2, 6]
    flipped = ring * np.array([1.0, -1.0, 1.0])            # mirror through the substituent axis
    rmsd = calculate_symmetry_corrected_rmsd(ring, flipped[perm], elements=el,
                                             elements_dock=[el[i] for i in perm])
    assert np.isclose(rmsd, 0.0, atol=1e-9)
    # the per-element Hungarian lower bound would also give 0 here; a genuine displacement must not
    shifted = ring + np.array([0.0, 0.0, 1.5])
    assert np.isclose(calculate_symmetry_corrected_rmsd(ring, shifted[perm], elements=el,
                                                        elements_dock=[el[i] for i in perm]), 1.5)
