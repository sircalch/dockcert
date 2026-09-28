"""
Rule-based receptor PDBQT writer shared by the benchmark scripts.

AutoDock Vina derives its own atom types (hydrophobic/polar carbon, donor, acceptor, metal) from the
AutoDock types in the PDBQT file and the bonding, so only the following decisions are needed:
  H   polar (within 1.35 A of N, O or S; PDBFixer places H at ~1.0-1.2 A) -> HD; non-polar H dropped
  O   -> OA (acceptor)
  S   -> SA
  N   -> NA if it carries no H and has exactly two heavy-atom neighbours (sp2 N with a lone pair:
         unprotonated His, adenine N1/N3/N7, ...); otherwise N
  C   -> C (Vina does not distinguish aromatic carbon)
  metals / halides: Zn, Fe, Mg, Mn, Ca, Cl, Br, I, F, P kept as such; elements without an AutoDock
  Vina 1.2 type (e.g. V, Na, K, Hg) are left out and reported
Partial charges are written as zero; Vina does not use them.
"""
import re

import numpy as np
from scipy.spatial import cKDTree

VINA_METALS = {"ZN": "Zn", "FE": "Fe", "MG": "Mg", "MN": "Mn", "CA": "Ca", "CL": "Cl", "BR": "Br", "I": "I",
               "F": "F", "P": "P"}
METALS = {"ZN", "FE", "MG", "MN", "CA"}
OTHER_METALS = {"NA", "K", "LI", "CO", "NI", "CU", "CD", "HG", "SR", "BA", "CS", "RB", "AL", "GA", "PT", "AU", "AG", "PB", "YB", "SM", "GD", "W", "MO", "V", "CR"}


def element_of(line):
    el = line[76:78].strip().upper() if len(line) >= 78 else ""
    if not el:
        name = line[12:16].strip()
        el = re.sub(r"[^A-Za-z]", "", name)[:1].upper()
    return el


def write_receptor_pdbqt(atom_lines, out_path):
    """atom_lines: PDB ATOM/HETATM records (hydrogens optional). Returns a count of the types written."""
    els = [element_of(l) for l in atom_lines]
    xyz = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])] for l in atom_lines])
    tree = cKDTree(xyz)
    heavy_nb = np.zeros(len(els), int)
    h_nb = np.zeros(len(els), int)
    polar_h = set()
    for i, j in tree.query_pairs(1.9):
        d = float(np.linalg.norm(xyz[i] - xyz[j]))
        for a, b in ((i, j), (j, i)):
            if els[a] == "H":
                if els[b] in ("N", "O", "S") and d < 1.35:
                    polar_h.add(a)
            elif els[b] == "H":
                if d < 1.35:
                    h_nb[a] += 1
            elif d < 1.9 and not ({els[a], els[b]} & (METALS | OTHER_METALS)):
                heavy_nb[a] += 1
    counts = {}
    out = []
    dropped = []
    for k, (line, e) in enumerate(zip(atom_lines, els)):
        if e == "H":
            if k not in polar_h:
                continue
            t = "HD"
        elif e == "O":
            t = "OA"
        elif e == "S":
            t = "SA"
        elif e == "N":
            t = "NA" if (h_nb[k] == 0 and heavy_nb[k] == 2) else "N"
        elif e == "C":
            t = "C"
        elif e in VINA_METALS:
            t = VINA_METALS[e]
        elif e in OTHER_METALS or e not in ("C", "N", "O", "S", "H"):
            dropped.append(line[17:20].strip() + ":" + e)   # no AutoDock Vina 1.2 type; left out
            continue
        else:
            t = e.capitalize()
        counts[t] = counts.get(t, 0) + 1
        out.append(f"{line[:54].ljust(54)}  1.00  0.00    +0.000 {t:<2}")
    open(out_path, "w").write("\n".join(out) + "\n")
    if dropped:
        counts["dropped"] = sorted(set(dropped))
    return counts
