"""
Tests for CSV and structure parsers in DockCert.
"""

import os
import tempfile
import numpy as np
import pytest
from dockcert.parsers.generic_csv import load_docking_csv
from dockcert.parsers.vina_smina import parse_vina_log
from dockcert.parsers.structure_io import load_molecule_coordinates


def test_generic_csv_parser():
    csv_content = """LigandID,docking_score,is_active,MW
CHEMBL1,-9.5,active,350.2
CHEMBL2,-8.1,active,360.1
CHEMBL3,-5.2,decoy,355.0
CHEMBL4,-4.8,decoy,348.0
"""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
        f.write(csv_content)
        f_path = f.name
        
    try:
        labels, scores, props, meta = load_docking_csv(f_path)
        assert len(labels) == 4
        assert labels[0] == 1
        assert labels[2] == 0
        assert scores[0] == -9.5
        assert "MW" in props
        assert meta["n_actives"] == 2
        assert meta["n_decoys"] == 2
    finally:
        if os.path.exists(f_path):
            os.remove(f_path)


def test_vina_log_parser():
    vina_log = """
#################################################################
# If you used AutoDock Vina in your work, please cite:          #
#################################################################

mode |   affinity | dist from best mode
     | (kcal/mol) | rmsd l.b.| rmsd u.b.
-----+------------+----------+----------
   1         -9.2      0.000      0.000
   2         -8.8      1.450      1.920
   3         -7.5      2.300      3.100
"""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".log", delete=False) as f:
        f.write(vina_log)
        f_path = f.name
        
    try:
        affs, lb, ub = parse_vina_log(f_path)
        assert len(affs) == 3
        assert affs[0] == -9.2
        assert lb[1] == 1.450
    finally:
        if os.path.exists(f_path):
            os.remove(f_path)


def test_pdb_coordinate_parser():
    pdb_content = """ATOM      1  C1  LIG A   1      10.000  20.000  30.000  1.00 20.00           C
ATOM      2  O1  LIG A   1      11.000  20.000  30.000  1.00 20.00           O
ATOM      3  H1  LIG A   1      10.000  21.000  30.000  1.00 20.00           H
"""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".pdb", delete=False) as f:
        f.write(pdb_content)
        f_path = f.name
        
    try:
        coords, elements = load_molecule_coordinates(f_path)
        assert len(coords) == 2  # Hydrogen ignored
        assert elements == ["C", "O"]
        assert coords[0, 0] == 10.0
    finally:
        if os.path.exists(f_path):
            os.remove(f_path)
