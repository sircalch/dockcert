"""
Real redocking with AutoDock Vina 1.2.7 (validation/redock_vina.py): trypsin-benzamidine (3PTB),
streptavidin-biotin (1STP) and HIV-1 protease-indinavir (1HSG). The expected RMSDs were computed
independently with RDKit rdMolAlign.CalcRMS on meeko-reconstructed poses.
"""
import os
import numpy as np
import pytest
from dockcert.parsers.structure_io import load_poses, load_molecule_coordinates
from dockcert.parsers.vina_smina import parse_vina_log
from dockcert.core.rmsd import calculate_symmetry_corrected_rmsd, evaluate_redocking_success

D = os.path.join(os.path.dirname(__file__), "data", "redock")
# RDKit CalcRMS (symmetry-aware, no superposition) for the 9 ranked poses
RDKIT = {
    "3PTB": [0.38, 0.41, 3.38, 3.62, 4.04, 11.16, 11.01, 3.47, 9.86],
    "1STP": [0.60, 1.97, 2.13, 2.37, 6.67, 8.10, 6.17, 6.80, 6.56],
    "1HSG": [4.39, 10.21, 10.85, 4.11, 10.80, 10.63, 10.10, 11.21, 10.77],
}
N_HEAVY = {"3PTB": 9, "1STP": 16, "1HSG": 45}


@pytest.mark.parametrize("pid", sorted(RDKIT))
def test_rmsd_matches_rdkit(pid):
    c_ref, el_ref = load_molecule_coordinates(os.path.join(D, f"{pid}_ref.sdf"))
    poses = load_poses(os.path.join(D, f"{pid}_out.pdbqt"))
    assert len(poses) == 9
    assert all(len(c) == N_HEAVY[pid] for c, _ in poses)          # hydrogens (HD) dropped
    assert all("A" not in el and "OA" not in el for _, el in poses)  # AutoDock types mapped to elements
    got = [calculate_symmetry_corrected_rmsd(c_ref, c, elements=el_ref, elements_dock=el) for c, el in poses]
    np.testing.assert_allclose(got, RDKIT[pid], atol=0.006)


def test_redocking_verdicts():
    status = {}
    for pid, vals in RDKIT.items():
        status[pid] = evaluate_redocking_success(vals, ranked_poses=True)["status"]
    assert status == {"3PTB": "PASS", "1STP": "PASS", "1HSG": "FAIL"}


def test_vina_log_real():
    aff, lb, ub = parse_vina_log(os.path.join(D, "3PTB_vina.log"))
    assert len(aff) == 9 and aff[0] == pytest.approx(-6.087) and lb[0] == 0.0


def test_macrocycle_pseudo_atoms_are_skipped():
    """meeko closes flexible rings with G0/CG0 pairs; G0 is a dummy atom, CG0 a real carbon."""
    d = os.path.join(os.path.dirname(__file__), "data", "redock_macrocycle")
    c_ref, el_ref = load_molecule_coordinates(os.path.join(d, "1MZC_BNE_ligand.sdf"))
    poses = load_poses(os.path.join(d, "1MZC_BNE_out.pdbqt"))
    assert all(len(c) == len(c_ref) == 35 for c, _ in poses)
    r = calculate_symmetry_corrected_rmsd(c_ref, poses[0][0], elements=el_ref, elements_dock=poses[0][1])
    assert np.isfinite(r)
