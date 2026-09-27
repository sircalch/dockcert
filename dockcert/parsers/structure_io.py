"""
Coordinate extraction from SDF, PDB, and PDBQT molecular structure files.

Multi-pose files are split into poses: MODEL/ENDMDL blocks in PDB/PDBQT (the Vina, Smina and
GNINA output format) and $$$$-separated records in SDF. In PDBQT files the element is derived
from the AutoDock atom type (A -> C, OA -> O, NA -> N, SA -> S, HD/HS -> H, ...), not read
verbatim from the type column.
"""

from typing import Tuple, List
import os
import numpy as np

# AutoDock 4 / Vina atom types -> element
_AD_TYPES = {
    "A": "C", "C": "C", "N": "N", "NA": "N", "NS": "N", "OA": "O", "OS": "O", "O": "O",
    "S": "S", "SA": "S", "H": "H", "HD": "H", "HS": "H", "P": "P", "F": "F", "CL": "CL",
    "Cl": "CL", "BR": "BR", "Br": "BR", "I": "I", "B": "B", "SE": "SE", "Se": "SE",
    "G0": "C", "G1": "C", "G2": "C", "G3": "C", "CG0": "C", "CG1": "C", "CG2": "C", "CG3": "C",
    "W": "O", "MG": "MG", "Mg": "MG", "ZN": "ZN", "Zn": "ZN", "CA": "CA", "Ca": "CA",
    "FE": "FE", "Fe": "FE", "MN": "MN", "Mn": "MN",
}

Pose = Tuple[np.ndarray, List[str]]


def _pdb_element(line: str, is_pdbqt: bool) -> str:
    if is_pdbqt:
        ad_type = line[77:79].strip() or line[76:79].strip()
        if ad_type in _AD_TYPES:
            return _AD_TYPES[ad_type]
        if ad_type.upper() in _AD_TYPES:
            return _AD_TYPES[ad_type.upper()]
    elem = line[76:78].strip()
    if not elem:
        name = line[12:16].strip()
        elem = "".join(c for c in name if c.isalpha())[:1]
    return elem.upper()


def _parse_pdb_poses(filepath: str, is_pdbqt: bool) -> List[Pose]:
    poses: List[Pose] = []
    coords: List[List[float]] = []
    elements: List[str] = []
    with open(filepath, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith("MODEL"):
                coords, elements = [], []
            elif line.startswith("ENDMDL"):
                if coords:
                    poses.append((np.array(coords, dtype=np.float64), elements))
                coords, elements = [], []
            elif line.startswith(("HETATM", "ATOM")):
                try:
                    xyz = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
                except ValueError:
                    continue
                elem = _pdb_element(line, is_pdbqt)
                if elem != "H":                                  # heavy atoms only
                    coords.append(xyz)
                    elements.append(elem)
    if coords:                                                   # file without MODEL records
        poses.append((np.array(coords, dtype=np.float64), elements))
    return poses


def _parse_sdf_poses(filepath: str) -> List[Pose]:
    with open(filepath, "r", encoding="utf-8", errors="ignore") as f:
        records = f.read().split("$$$$")
    poses: List[Pose] = []
    for k, rec in enumerate(records):
        lines = rec.splitlines()
        if k > 0 and rec.startswith(("\r", "\n")):          # rest of the "$$$$" line (the title line may be blank)
            lines = lines[1:]
        if len(lines) < 4 or not "".join(lines).strip():
            continue
        try:
            n_atoms = int(lines[3][:3])
        except ValueError:
            continue
        coords, elements = [], []
        for ln in lines[4:4 + n_atoms]:
            elem = ln[31:34].strip().upper()
            if elem != "H":
                coords.append([float(ln[:10]), float(ln[10:20]), float(ln[20:30])])
                elements.append(elem)
        if coords:
            poses.append((np.array(coords, dtype=np.float64), elements))
    return poses


def load_poses(filepath: str) -> List[Pose]:
    """
    Returns every pose in a structure file as (heavy-atom coordinates [N, 3], element symbols).

    For docking output the poses keep the program's ranking order (pose 0 = top-scored).
    """
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"Structure file not found: {filepath}")
    ext = os.path.splitext(filepath)[1].lower()
    if ext in (".pdb", ".pdbqt", ".ent"):
        poses = _parse_pdb_poses(filepath, is_pdbqt=(ext == ".pdbqt"))
    elif ext in (".sdf", ".mol"):
        poses = _parse_sdf_poses(filepath)
    else:
        raise ValueError(f"Unsupported structure format: {ext}")
    if not poses:
        raise ValueError(f"Could not extract 3D coordinates from {filepath}")
    return poses


def load_molecule_coordinates(filepath: str) -> Tuple[np.ndarray, List[str]]:
    """
    Extracts 3D heavy-atom coordinates and element symbols of the FIRST pose of an SDF, PDB,
    or PDBQT file. Use load_poses() for multi-pose docking output.
    """
    return load_poses(filepath)[0]
