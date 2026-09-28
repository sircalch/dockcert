"""
Redocking benchmark on the Astex Diverse Set (85 complexes; Hartshorn et al., J. Med. Chem. 2007)
as distributed with PoseBusters (Buttenschoen et al., Chem. Sci. 2024; Zenodo 10.5281/zenodo.8278563,
CC-BY 4.0), with AutoDock Vina 1.2.7.

For every complex:
  receptor  <id>_protein.pdb -> PDBFixer (modified residues, missing atoms, hydrogens at pH 7.4);
            waters removed; cofactors, metals and ions within 8 A of the ligand kept (heavy atoms);
            typed by pdbqt_rules.py
  ligand    <id>_ligand_start_conf.sdf (RDKit conformer supplied by PoseBusters, not the crystal pose)
            -> meeko MoleculePreparation (PDBQT)
  box       centred on the crystal ligand, edge = max(22 A, ligand extent + 10 A)
  Vina      exhaustiveness 8 (default), 9 modes, seed 42

Usage:  python astex_redock.py DATA_DIR OUT_DIR [--workers 4] [--cpu 2] [--only ID ...]
DATA_DIR is the astex_diverse_set folder from posebusters_paper_data.zip.
Results: OUT_DIR/<id>/{rec.pdbqt, lig.pdbqt, out.pdbqt, vina.log, status.json}
"""
import argparse
import json
import os
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np

VINA_DEFAULT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "vina.exe")


WATER = {"HOH", "WAT", "DOD"}


def prep_receptor(pdb_in, out_dir, ligand_sdf, het_cutoff=8.0):
    """
    Protein: PDBFixer (modified residues -> standard, missing heavy atoms, hydrogens at pH 7.4).
    Heterogens: non-water HETATM residues with any atom within het_cutoff A of the crystal ligand
    (cofactors, metals, ions) are kept as heavy atoms. Everything is typed by pdbqt_rules.
    """
    from pdbfixer import PDBFixer
    from openmm.app import PDBFile
    from rdkit import Chem
    from pdbqt_rules import write_receptor_pdbqt
    fx = PDBFixer(filename=pdb_in)
    fx.findNonstandardResidues()
    fx.replaceNonstandardResidues()
    fx.removeHeterogens(keepWater=False)
    fx.findMissingResidues()
    fx.missingResidues = {}
    fx.findMissingAtoms()
    fx.addMissingAtoms()
    fx.addMissingHydrogens(7.4)
    h_pdb = os.path.join(out_dir, "recH.pdb")
    with open(h_pdb, "w") as fh:
        PDBFile.writeFile(fx.topology, fx.positions, fh, keepIds=True)
    protein = [l for l in open(h_pdb) if l.startswith(("ATOM", "HETATM"))]

    lig = Chem.MolFromMolFile(ligand_sdf, removeHs=True).GetConformer().GetPositions()
    standard = {l[17:20] for l in protein}
    het_by_res = {}
    for l in open(pdb_in):
        if l.startswith("HETATM") and l[17:20].strip() not in WATER and l[17:20] not in standard:
            if l[76:78].strip().upper() == "H":
                continue
            het_by_res.setdefault((l[17:20], l[21], l[22:27]), []).append(l)
    kept = []
    for key, lines in het_by_res.items():
        xyz = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])] for l in lines])
        if np.min(np.linalg.norm(xyz[:, None, :] - lig[None, :, :], axis=-1)) <= het_cutoff:
            kept.append(key[0].strip())
            protein += lines
    counts = write_receptor_pdbqt(protein, os.path.join(out_dir, "rec.pdbqt"))
    return sorted(kept), counts


def prep_ligand(sdf_in, out_dir):
    from rdkit import Chem
    from meeko import MoleculePreparation, PDBQTWriterLegacy
    mol = Chem.MolFromMolFile(sdf_in, removeHs=False)
    if mol.GetNumAtoms() == Chem.RemoveHs(mol).GetNumAtoms():
        mol = Chem.AddHs(mol, addCoords=True)
    setups = MoleculePreparation().prepare(mol)
    pdbqt, ok, err = PDBQTWriterLegacy.write_string(setups[0])
    if not ok:
        raise RuntimeError("ligand preparation failed: " + str(err))
    open(os.path.join(out_dir, "lig.pdbqt"), "w").write(pdbqt)


def box_from_crystal(sdf_crystal):
    from rdkit import Chem
    m = Chem.MolFromMolFile(sdf_crystal, removeHs=True)
    xyz = m.GetConformer().GetPositions()
    centre = xyz.mean(axis=0)
    size = max(22.0, float(np.ptp(xyz, axis=0).max()) + 10.0)
    return centre, size


def run_one(args):
    cid, data_dir, out_root, vina, cpu = args
    src = os.path.join(data_dir, cid)
    out = os.path.join(out_root, cid)
    os.makedirs(out, exist_ok=True)
    status_p = os.path.join(out, "status.json")
    if os.path.exists(status_p) and json.load(open(status_p)).get("ok"):
        return cid, "cached"
    t0 = time.time()
    st = {"id": cid, "ok": False}
    try:
        kept, counts = prep_receptor(os.path.join(src, f"{cid}_protein.pdb"), out,
                                     os.path.join(src, f"{cid}_ligand.sdf"))
        st["kept_heterogens"] = kept
        if "dropped" in counts:
            st["dropped_atoms"] = counts["dropped"]
        prep_ligand(os.path.join(src, f"{cid}_ligand_start_conf.sdf"), out)
        c, size = box_from_crystal(os.path.join(src, f"{cid}_ligand.sdf"))
        cmd = [vina, "--receptor", "rec.pdbqt", "--ligand", "lig.pdbqt",
               "--center_x", f"{c[0]:.3f}", "--center_y", f"{c[1]:.3f}", "--center_z", f"{c[2]:.3f}",
               "--size_x", f"{size:.1f}", "--size_y", f"{size:.1f}", "--size_z", f"{size:.1f}",
               "--exhaustiveness", "8", "--num_modes", "9", "--seed", "42", "--cpu", str(cpu),
               "--out", "out.pdbqt"]
        r = subprocess.run(cmd, cwd=out, capture_output=True, text=True)
        open(os.path.join(out, "vina.log"), "w").write(r.stdout + r.stderr)
        if r.returncode != 0 or not os.path.exists(os.path.join(out, "out.pdbqt")):
            raise RuntimeError("vina failed: " + (r.stderr or r.stdout)[-400:])
        st.update(ok=True, box_centre=[float(x) for x in c], box_size=size, seconds=round(time.time() - t0, 1))
    except Exception as e:                                   # recorded, not hidden
        st["error"] = str(e)[:600]
    json.dump(st, open(status_p, "w"), indent=1)
    return cid, "ok" if st["ok"] else "FAILED: " + st.get("error", "")[:120]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("data_dir")
    ap.add_argument("out_dir")
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--cpu", type=int, default=2)
    ap.add_argument("--vina", default=VINA_DEFAULT)
    ap.add_argument("--only", nargs="*")
    a = ap.parse_args()
    ids = sorted(d for d in os.listdir(a.data_dir) if os.path.isdir(os.path.join(a.data_dir, d)))
    if a.only:
        ids = [i for i in ids if i in a.only]
    jobs = [(i, a.data_dir, a.out_dir, a.vina, a.cpu) for i in ids]
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        for cid, msg in ex.map(run_one, jobs):
            print(cid, msg, flush=True)


if __name__ == "__main__":
    main()
