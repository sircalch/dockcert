"""
Coordinate extraction from SDF, PDB, and PDBQT molecular structure files.
"""

from typing import Tuple, List, Optional
import os
import numpy as np


def load_molecule_coordinates(filepath: str) -> Tuple[np.ndarray, List[str]]:
    """
    Extracts 3D heavy-atom coordinates and element symbols from SDF, PDB, or PDBQT.

    Parameters
    ----------
    filepath : str
        Path to structural coordinate file.

    Returns
    -------
    coords : np.ndarray
        Array of (x, y, z) coordinates (shape: [N_atoms, 3]).
    elements : list of str
        Element symbol for each atom.
    """
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"Structure file not found: {filepath}")
        
    ext = os.path.splitext(filepath)[1].lower()
    
    coords = []
    elements = []
    
    if ext in [".pdb", ".pdbqt", ".ent"]:
        with open(filepath, "r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                if line.startswith("HETATM") or line.startswith("ATOM"):
                    try:
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                        # Element symbol usually in columns 76-78 or derived from atom name (cols 12-16)
                        elem = line[76:78].strip()
                        if not elem:
                            name = line[12:16].strip()
                            elem = "".join([c for c in name if c.isalpha()])[:2].capitalize()
                            
                        # Ignore hydrogens for standard heavy-atom RMSD
                        if elem.upper() != "H":
                            coords.append([x, y, z])
                            elements.append(elem.upper())
                    except (ValueError, IndexError):
                        continue
                        
    elif ext in [".sdf", ".mol"]:
        with open(filepath, "r", encoding="utf-8", errors="ignore") as f:
            lines = f.readlines()
            
        if len(lines) >= 4:
            # Counts line is line index 3 (4th line)
            counts_line = lines[3]
            try:
                n_atoms = int(counts_line[:3])
                for i in range(4, 4 + n_atoms):
                    atom_line = lines[i]
                    x = float(atom_line[:10])
                    y = float(atom_line[10:20])
                    z = float(atom_line[20:30])
                    elem = atom_line[31:34].strip().upper()
                    
                    if elem != "H":
                        coords.append([x, y, z])
                        elements.append(elem)
            except (ValueError, IndexError):
                pass

    if not coords:
        raise ValueError(f"Could not extract 3D coordinates from {filepath}")
        
    return np.array(coords, dtype=np.float64), elements
