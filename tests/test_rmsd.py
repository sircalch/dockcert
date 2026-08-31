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


def test_evaluate_redocking_success():
    res_pass = evaluate_redocking_success([1.2, 1.4, 2.5, 3.2])
    assert res_pass["status"] == "PASS"
    assert res_pass["min_rmsd"] == 1.2
    
    res_fail = evaluate_redocking_success([4.5, 5.2, 6.1])
    assert res_fail["status"] == "FAIL"
