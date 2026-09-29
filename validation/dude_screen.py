"""
Retrospective virtual screen of DUD-E targets (Mysinger et al., J. Med. Chem. 2012) with AutoDock
Vina 1.2.7, to exercise DockCert's enrichment statistics on real docking scores.

Receptor  DUD-E receptor.pdb (polar hydrogens present, DOCK residue names, no element column; the
          element is taken from the first letter of the atom name) -> PDBQT with the rules of
          pdbqt_rules.py (same writer as the Astex benchmark).
Ligands   actives_final.ism / decoys_final.ism (protonation as supplied by DUD-E) -> RDKit
          ETKDGv3 + MMFF94 conformer (seed 42) -> meeko PDBQT
Box       centred on crystal_ligand.mol2, edge = max(22 A, ligand extent + 10 A)
Vina      exhaustiveness 8, 1 mode, seed 42, one CPU per ligand

Usage:  python dude_screen.py TARGET_DIR OUT_CSV [--workers 8]
OUT_CSV is appended ligand by ligand (restartable); failures are recorded with their reason.
"""
import argparse
import csv
import os
import re
import subprocess
import tempfile
from concurrent.futures import ProcessPoolExecutor

import numpy as np

VINA_DEFAULT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "vina.exe")


def receptor_pdbqt(pdb_path, out_path):
    from pdbqt_rules import write_receptor_pdbqt
    lines = [l.rstrip("\r\n") for l in open(pdb_path) if l.startswith(("ATOM", "HETATM"))]
    return write_receptor_pdbqt(lines, out_path)


def box_from_mol2(path):
    xyz, on = [], False
    for line in open(path):
        if line.startswith("@<TRIPOS>ATOM"):
            on = True
            continue
        if line.startswith("@<TRIPOS>") and on:
            break
        if on and line.strip():
            p = line.split()
            if not p[5].upper().startswith("H"):
                xyz.append([float(p[2]), float(p[3]), float(p[4])])
    xyz = np.array(xyz)
    return xyz.mean(axis=0), max(22.0, float(np.ptp(xyz, axis=0).max()) + 10.0)


def ligand_pdbqt(smiles):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from meeko import MoleculePreparation, PDBQTWriterLegacy
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        raise ValueError("SMILES not parsed")
    m = Chem.AddHs(m)
    p = AllChem.ETKDGv3()
    p.randomSeed = 42
    if AllChem.EmbedMolecule(m, p) != 0:
        raise ValueError("embedding failed")
    AllChem.MMFFOptimizeMolecule(m, maxIters=2000)
    setup = MoleculePreparation().prepare(m)[0]
    s, ok, err = PDBQTWriterLegacy.write_string(setup)
    if not ok:
        raise ValueError(f"meeko: {err}")
    return s


def dock(job):
    lid, label, smiles, rec, centre, size, vina = job
    try:
        lig = ligand_pdbqt(smiles)
        with tempfile.TemporaryDirectory() as td:
            lp = os.path.join(td, "lig.pdbqt")
            open(lp, "w").write(lig)
            cmd = [vina, "--receptor", rec, "--ligand", lp,
                   "--center_x", f"{centre[0]:.3f}", "--center_y", f"{centre[1]:.3f}", "--center_z", f"{centre[2]:.3f}",
                   "--size_x", f"{size:.1f}", "--size_y", f"{size:.1f}", "--size_z", f"{size:.1f}",
                   "--exhaustiveness", "8", "--num_modes", "1", "--seed", "42", "--cpu", "1",
                   "--out", os.path.join(td, "out.pdbqt")]
            for _attempt in range(2):                 # one retry for transient failures
                r = subprocess.run(cmd, capture_output=True, text=True)
                m = re.search(r"^\s+1\s+(-?\d+(?:\.\d+)?)\s", r.stdout, re.M)   # Vina writes -7, not -7.000
                if m:
                    break
            if not m:
                raise RuntimeError("no score: " + (r.stderr or r.stdout)[-200:])
            return lid, label, float(m.group(1)), ""
    except Exception as e:
        return lid, label, "", str(e)[:150].replace("\n", " ")


def keep_awake():
    """Ask Windows not to sleep while this process runs (released automatically when it exits)."""
    if os.name == "nt":
        import ctypes
        ES_CONTINUOUS, ES_SYSTEM_REQUIRED = 0x80000000, 0x00000001
        ctypes.windll.kernel32.SetThreadExecutionState(ES_CONTINUOUS | ES_SYSTEM_REQUIRED)


def main():
    keep_awake()
    ap = argparse.ArgumentParser()
    ap.add_argument("target_dir")
    ap.add_argument("out_csv")
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--vina", default=VINA_DEFAULT)
    a = ap.parse_args()
    rec = os.path.join(a.target_dir, "receptor.pdbqt")
    receptor_pdbqt(os.path.join(a.target_dir, "receptor.pdb"), rec)
    centre, size = box_from_mol2(os.path.join(a.target_dir, "crystal_ligand.mol2"))
    ligs = []
    for fname, label in (("actives_final.ism", 1), ("decoys_final.ism", 0)):
        for k, line in enumerate(open(os.path.join(a.target_dir, fname))):
            p = line.split()
            if p:
                ligs.append((f"{'A' if label else 'D'}{k:05d}_{p[1] if len(p) > 1 else ''}", label, p[0]))
    done = set()
    if os.path.exists(a.out_csv):
        # rows with a score or with a permanent error count as done; transient Vina failures
        # (no score line, e.g. a run interrupted by system sleep) are retried
        done = {r["id"] for r in csv.DictReader(open(a.out_csv))
                if r["vina_score"] or not r["error"].startswith("no score")}
    new = not os.path.exists(a.out_csv)
    jobs = [(i, lab, smi, rec, centre, size, a.vina) for i, lab, smi in ligs if i not in done]
    with open(a.out_csv, "a", newline="") as fh:
        w = csv.writer(fh)
        if new:
            w.writerow(["id", "label", "vina_score", "error"])
        with ProcessPoolExecutor(max_workers=a.workers) as ex:
            for n, row in enumerate(ex.map(dock, jobs, chunksize=4)):
                w.writerow(row)
                if n % 50 == 0:
                    fh.flush()
    print("done", a.target_dir)


if __name__ == "__main__":
    main()
