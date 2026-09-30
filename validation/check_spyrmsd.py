"""
Second independent reference for the symmetry-corrected RMSD: spyrmsd (Meli & Biggin 2020), which
matches atoms by graph isomorphism on connectivity and element only, as DockCert does.

    python validation/check_spyrmsd.py DATA_DIR RUN_DIR

DATA_DIR and RUN_DIR are those of analyze_astex.py (PoseBusters Astex files; Vina runs, as in
results/astex_vina_outputs.zip). Adds rmsd_spyrmsd to results/astex_poses.csv and prints the agreement.
"""
import os
import sys

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from meeko import PDBQTMolecule, RDKitMolCreate
from spyrmsd import rmsd as sp_rmsd
from spyrmsd import molecule as sp_mol

RDLogger.DisableLog("rdApp.*")
HERE = os.path.dirname(os.path.abspath(__file__))


def spyrmsd_case(cid, data_dir, run_dir):
    crystal = Chem.MolFromMolFile(os.path.join(data_dir, cid, f"{cid}_ligand.sdf"), removeHs=True, sanitize=False)
    topo = Chem.RemoveAllHs(Chem.MolFromMolFile(os.path.join(data_dir, cid, f"{cid}_ligand_start_conf.sdf"), removeHs=False))
    ref = Chem.Mol(topo)
    ref.RemoveAllConformers()
    ref.AddConformer(Chem.Conformer(crystal.GetConformer()), assignId=True)
    poses = Chem.RemoveAllHs(RDKitMolCreate.from_pdbqt_mol(
        PDBQTMolecule.from_file(os.path.join(run_dir, cid, "out.pdbqt"), skip_typing=True))[0])
    mref = sp_mol.Molecule.from_rdkit(ref)
    mpose = sp_mol.Molecule.from_rdkit(poses)
    coords = [c.GetPositions() for c in poses.GetConformers()]
    # same topology for reference and poses; no centring or superposition (docking RMSD)
    return sp_rmsd.symmrmsd(mref.coordinates, coords, mref.atomicnums, mpose.atomicnums,
                            mref.adjacency_matrix, mpose.adjacency_matrix, center=False, minimize=False)


def main(data_dir, run_dir):
    p = pd.read_csv(os.path.join(HERE, "results", "astex_poses.csv"))
    vals = {}
    for cid in sorted(p.id.unique()):
        r = spyrmsd_case(cid, data_dir, run_dir)
        for k, v in enumerate(np.atleast_1d(r)):
            vals[(cid, k + 1)] = float(v)
    p["rmsd_spyrmsd"] = [vals.get((i, k), np.nan) for i, k in zip(p.id, p["rank"])]
    p.to_csv(os.path.join(HERE, "results", "astex_poses.csv"), index=False)
    d = (p.rmsd_dockcert - p.rmsd_spyrmsd).abs()
    print(f"{p.rmsd_spyrmsd.notna().sum()} poses; max |DockCert - spyrmsd| = {d.max():.2e} A; "
          f"{int((d <= 2e-6).sum())} within 2e-6 A")
    print("1TT1:", p[p.id.str.startswith("1TT1")][["rmsd_rdkit", "rmsd_dockcert", "rmsd_spyrmsd"]].round(4).to_string())


if __name__ == "__main__":
    main(*sys.argv[1:3])
