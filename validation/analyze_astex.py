"""
Pose-accuracy analysis of the Astex redocking campaign (astex_redock.py).

For every complex and every Vina pose (rank order) the heavy-atom RMSD to the crystal ligand, without
superposition, is computed in four ways:
  rdkit     RDKit rdMolAlign.CalcRMS on meeko-reconstructed poses: symmetry-aware graph matching
            (reference implementation)
  dockcert  DockCert 1.1: calculate_symmetry_corrected_rmsd on load_poses() (graph isomorphisms)
  naive     RMSD for one fixed atom correspondence (the first substructure match), i.e. without
            symmetry correction
  hung_lb   per-element Hungarian assignment (bonding ignored): a lower bound of the RMSD
  hung_v100 the Hungarian routine of DockCert 1.0.0 as it was called, with the reference element
            list applied to the pose atoms in their own order
and the single number that the DockCert 1.0.0 command line produced for the whole output file
(v100_cli: all poses read as one molecule, AutoDock types taken as elements).

Usage: python analyze_astex.py DATA_DIR RUN_DIR OUT_DIR
"""
import json
import os
import sys

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolAlign
from meeko import PDBQTMolecule, RDKitMolCreate

RDLogger.DisableLog("rdApp.*")
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)
from dockcert.parsers.structure_io import load_molecule_coordinates, load_poses          # noqa: E402
from dockcert.parsers.vina_smina import parse_vina_log                                    # noqa: E402
from dockcert.core.rmsd import calculate_symmetry_corrected_rmsd                          # noqa: E402
from legacy_v100 import rmsd_v100, structure_io_v100                                      # noqa: E402


def hungarian_lower_bound(c_ref, el_ref, c_pose, el_pose):
    from scipy.optimize import linear_sum_assignment
    tot = 0.0
    for e in set(el_ref):
        a = np.asarray(c_ref)[[i for i, x in enumerate(el_ref) if x == e]]
        b = np.asarray(c_pose)[[i for i, x in enumerate(el_pose) if x == e]]
        if len(a) != len(b):
            return np.nan
        cost = np.sum((a[:, None, :] - b[None, :, :]) ** 2, axis=-1)
        r, c = linear_sum_assignment(cost)
        tot += float(cost[r, c].sum())
    return float(np.sqrt(tot / len(el_ref)))


def analyse(cid, data_dir, run_dir):
    run = os.path.join(run_dir, cid)
    st = json.load(open(os.path.join(run, "status.json")))
    if not st.get("ok"):
        return [], {"id": cid, "ok": False, "error": st.get("error", "")}
    ref_sdf = os.path.join(data_dir, cid, f"{cid}_ligand.sdf")
    out = os.path.join(run, "out.pdbqt")
    # Reference: crystal coordinates on the topology of the supplied start conformer, so that bond
    # orders and charges agree with the docked molecule
    crystal = Chem.MolFromMolFile(ref_sdf, removeHs=True, sanitize=False)
    topo = Chem.RemoveAllHs(Chem.MolFromMolFile(os.path.join(data_dir, cid, f"{cid}_ligand_start_conf.sdf"), removeHs=False))
    if crystal.GetNumAtoms() != topo.GetNumAtoms():
        return [], {"id": cid, "ok": False, "error": f"heavy-atom count differs: crystal {crystal.GetNumAtoms()}, start conformer {topo.GetNumAtoms()}"}
    ref = Chem.Mol(topo)
    ref.RemoveAllConformers()
    ref.AddConformer(Chem.Conformer(crystal.GetConformer()), assignId=True)
    if [a.GetSymbol() for a in topo.GetAtoms()] != [a.GetSymbol() for a in crystal.GetAtoms()]:
        return [], {"id": cid, "ok": False, "error": "crystal and start-conformer atom orders differ"}
    affs = parse_vina_log(os.path.join(run, "vina.log"))[0]

    pm = PDBQTMolecule.from_file(out, skip_typing=True)
    mols = RDKitMolCreate.from_pdbqt_mol(pm)
    poses_rd = Chem.RemoveAllHs(mols[0])
    c_ref, el_ref = load_molecule_coordinates(ref_sdf)
    poses = load_poses(out)
    rows = []
    for k, conf in enumerate(poses_rd.GetConformers()):
        pose = Chem.Mol(poses_rd, confId=conf.GetId())
        pose.RemoveAllConformers()
        pose.AddConformer(Chem.Conformer(conf), assignId=True)
        r_rdkit = rdMolAlign.CalcRMS(pose, ref)
        match = pose.GetSubstructMatch(ref)                     # one fixed correspondence ref -> pose
        if match:
            xyz = conf.GetPositions()[list(match)]
            r_naive = float(np.sqrt(np.mean(np.sum((xyz - ref.GetConformer().GetPositions()) ** 2, axis=1))))
        else:
            r_naive = np.nan
        c, el = poses[k]
        r_dc = calculate_symmetry_corrected_rmsd(c_ref, c, elements=el_ref, elements_dock=el)
        try:
            r_hung = rmsd_v100.calculate_symmetry_corrected_rmsd(c_ref, c, elements=el_ref)
        except Exception:
            r_hung = np.nan
        r_lb = hungarian_lower_bound(c_ref, el_ref, c, el)
        rows.append({"id": cid, "rank": k + 1, "vina_kcal": affs[k] if k < len(affs) else np.nan,
                     "rmsd_rdkit": r_rdkit, "rmsd_dockcert": r_dc, "rmsd_naive": r_naive,
                     "rmsd_hung_v100": r_hung, "rmsd_hung_lb": r_lb, "n_heavy": ref.GetNumAtoms()})
    # what the 1.0.0 command line reported for this output file
    try:
        c1, e1 = structure_io_v100.load_molecule_coordinates(ref_sdf)
        c2, e2 = structure_io_v100.load_molecule_coordinates(out)
        v100 = rmsd_v100.calculate_symmetry_corrected_rmsd(c1, c2, elements=e1)
        v100_status = rmsd_v100.evaluate_redocking_success([v100])["status"]
    except Exception as e:
        v100, v100_status = np.nan, f"error: {e}"[:60]
    return rows, {"id": cid, "ok": True, "n_poses": len(rows), "v100_cli_rmsd": v100, "v100_cli_status": v100_status,
                  "box_size": st.get("box_size"), "seconds": st.get("seconds")}


def main(data_dir, run_dir, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    ids = sorted(d for d in os.listdir(run_dir) if os.path.isdir(os.path.join(run_dir, d)))
    all_rows, cases = [], []
    for cid in ids:
        try:
            rows, case = analyse(cid, data_dir, run_dir)
        except Exception as e:
            rows, case = [], {"id": cid, "ok": False, "error": f"analysis: {e}"[:200]}
        all_rows += rows
        cases.append(case)
    poses = pd.DataFrame(all_rows)
    cases = pd.DataFrame(cases)
    poses.to_csv(os.path.join(out_dir, "astex_poses.csv"), index=False)
    cases.to_csv(os.path.join(out_dir, "astex_cases.csv"), index=False)
    print(f"{int(cases.ok.sum())}/{len(cases)} complexes analysed; {len(poses)} poses")
    return poses, cases


if __name__ == "__main__":
    main(*sys.argv[1:4])
