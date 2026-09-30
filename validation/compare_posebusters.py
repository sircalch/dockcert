"""
Compares the Astex redocking of this benchmark with the published PoseBusters results for Vina.

    python validation/compare_posebusters.py

Downloads posebusters_paper_results.csv from the PoseBusters data deposit (Zenodo 8278563,
CC-BY 4.0; Buttenschoen et al. 2024) into validation/data/ and prints the top-ranked success of
Vina reported there, ours (RDKit CalcRMS on the same criterion) and the number of complexes with
the same outcome (RMSD <= 2 A or not).
"""
import os
import urllib.request

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
URL = "https://zenodo.org/api/records/8278563/files/posebusters_paper_results.csv/content"
DATA = os.path.join(HERE, "data", "posebusters_paper_results.csv")


def main():
    if not os.path.exists(DATA):
        os.makedirs(os.path.dirname(DATA), exist_ok=True)
        urllib.request.urlretrieve(URL, DATA)
    pb = pd.read_csv(DATA)
    pb = pb[(pb.dataset == "astex") & (pb["post-processing"] == "none") & (pb.method == "vina")][["pdb_id", "rmsd"]]
    ours = pd.read_csv(os.path.join(HERE, "results", "astex_poses.csv"))
    ours = ours[ours["rank"] == 1].assign(pdb_id=lambda d: d.id.str[:4])
    m = pb.merge(ours, on="pdb_id")
    same = ((m.rmsd <= 2) == (m.rmsd_rdkit <= 2)).sum()
    print(f"complexes: {len(m)}; PoseBusters Vina success: {(m.rmsd <= 2).sum()}; "
          f"this benchmark: {(m.rmsd_rdkit <= 2).sum()}; same outcome: {same}")


if __name__ == "__main__":
    main()
