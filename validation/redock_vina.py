"""Redocking benchmark for DockCert: crystal complexes -> Vina 1.2.7 -> PDBQT poses + reference ligands.

Needs rdkit, meeko 0.8, pdbfixer/openmm and the Vina 1.2.7 executable saved as vina.exe next to this
script; the RCSB entries 3PTB, 1STP and 1HSG are expected here as <ID>.pdb. Outputs are copied to
tests/data/redock/. Usage: python redock_vina.py [PDB IDs]
"""
import os
import subprocess
import sys

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from pdbfixer import PDBFixer
from openmm.app import PDBFile
from meeko import MoleculePreparation, PDBQTWriterLegacy

HERE = os.path.dirname(os.path.abspath(__file__))
VINA = os.path.join(HERE, "vina.exe")
DROP = {"1STP": {"133"}}
CASES = {
    "3PTB": ("BEN", "NC(=N)c1ccccc1"),
    "1STP": ("BTN", "OC(=O)CCCC[C@@H]1SC[C@@H]2NC(=O)N[C@H]12"),
    "1HSG": ("MK1", "CC(C)(C)NC(=O)[C@@H]1CN(Cc2cccnc2)CCN1C[C@@H](O)C[C@@H](Cc1ccccc1)C(=O)N[C@@H]1c2ccccc2C[C@@H]1O"),
}


def run(cmd):
    r = subprocess.run(cmd, cwd=HERE, capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stdout[-2000:], r.stderr[-2000:])
        raise SystemExit(f"failed: {cmd[0]}")
    return r.stdout


def prep(pdb_id, resn, smiles):
    lines = open(os.path.join(HERE, pdb_id + ".pdb")).read().splitlines()
    lig = [l for l in lines if l.startswith("HETATM") and l[17:20] == resn]
    first = (lig[0][21], lig[0][22:26])                      # one copy of the ligand
    lig = [l for l in lig if (l[21], l[22:26]) == first and l[16] in " A"]
    open(os.path.join(HERE, f"{pdb_id}_ref.pdb"), "w").write("\n".join(lig) + "\nEND\n")

    ref = Chem.MolFromPDBFile(os.path.join(HERE, f"{pdb_id}_ref.pdb"), removeHs=True)
    ref = AllChem.AssignBondOrdersFromTemplate(Chem.MolFromSmiles(smiles), ref)
    Chem.MolToMolFile(ref, os.path.join(HERE, f"{pdb_id}_ref.sdf"))
    center = ref.GetConformer().GetPositions().mean(axis=0)

    # ligand input: RDKit conformer from SMILES (not the crystal pose), so the test is honest
    m = Chem.AddHs(Chem.MolFromSmiles(smiles))
    AllChem.EmbedMolecule(m, randomSeed=7)
    AllChem.MMFFOptimizeMolecule(m)
    setup = MoleculePreparation().prepare(m)[0]
    pdbqt, ok, err = PDBQTWriterLegacy.write_string(setup)
    open(os.path.join(HERE, f"{pdb_id}_lig.pdbqt"), "w").write(pdbqt)

    fx = PDBFixer(filename=os.path.join(HERE, pdb_id + ".pdb"))
    fx.removeHeterogens(keepWater=False)
    fx.findMissingResidues(); fx.missingResidues = {}
    fx.findMissingAtoms(); fx.addMissingAtoms(); fx.addMissingHydrogens(7.0)
    with open(os.path.join(HERE, f"{pdb_id}_recH.pdb"), "w") as fh:
        PDBFile.writeFile(fx.topology, fx.positions, fh, keepIds=True)
    if pdb_id in DROP:   # terminal residues meeko cannot template, far from the site
        p = os.path.join(HERE, f"{pdb_id}_recH.pdb")
        keep = [l for l in open(p) if not (l.startswith(("ATOM", "HETATM")) and l[22:26].strip() in DROP[pdb_id])]
        open(p, "w").writelines(keep)
    run(["mk_prepare_receptor.exe", "--read_pdb", f"{pdb_id}_recH.pdb", "-o", f"{pdb_id}_rec",
         "-p", "--allow_bad_res"])
    return center


def dock(pdb_id, center):
    log = run([VINA, "--receptor", f"{pdb_id}_rec.pdbqt", "--ligand", f"{pdb_id}_lig.pdbqt",
               "--center_x", f"{center[0]:.3f}", "--center_y", f"{center[1]:.3f}", "--center_z", f"{center[2]:.3f}",
               "--size_x", "22", "--size_y", "22", "--size_z", "22", "--exhaustiveness", "16",
               "--num_modes", "9", "--seed", "42", "--cpu", "8", "--out", f"{pdb_id}_out.pdbqt"])
    open(os.path.join(HERE, f"{pdb_id}_vina.log"), "w").write(log)


if __name__ == "__main__":
    for pid, (resn, smi) in CASES.items():
        if len(sys.argv) > 1 and pid not in sys.argv[1:]:
            continue
        c = prep(pid, resn, smi)
        dock(pid, c)
        print(pid, "docked", flush=True)
