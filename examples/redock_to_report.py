"""
Worked example: real redocking with AutoDock Vina, and a screening validation.

Uses the files bundled with the repository:
  tests/data/redock/<PDB>_ref.sdf     crystallographic ligand
  tests/data/redock/<PDB>_out.pdbqt   nine ranked poses from AutoDock Vina 1.2.7
  benchmark_screening.csv             a screening table (1,550 compounds, actives and decoys)

Run from the repository root:

    python examples/redock_to_report.py
"""

import os

from dockcert import assess_docking_quality
from dockcert.core.rmsd import calculate_symmetry_corrected_rmsd, evaluate_redocking_success
from dockcert.parsers.generic_csv import load_docking_csv
from dockcert.parsers.structure_io import load_molecule_coordinates, load_poses
from dockcert.reporters import (
    generate_docking_figures,
    generate_docking_html_report,
    generate_docking_manuscript_assets,
)

REDOCK_DIR = os.path.join("tests", "data", "redock")
SCREENING_CSV = "benchmark_screening.csv"
OUTPUT_DIR = "redock_example_output"
COMPLEXES = {"3PTB": "trypsin-benzamidine", "1STP": "streptavidin-biotin", "1HSG": "HIV-1 protease-indinavir"}


def redock_rmsds(pdb_id):
    """RMSD of every ranked pose to the crystallographic ligand (symmetry-aware, no superposition)."""
    ref_xyz, ref_elements = load_molecule_coordinates(os.path.join(REDOCK_DIR, f"{pdb_id}_ref.sdf"))
    poses = load_poses(os.path.join(REDOCK_DIR, f"{pdb_id}_out.pdbqt"))
    return [
        calculate_symmetry_corrected_rmsd(ref_xyz, xyz, elements=ref_elements, elements_dock=elements)
        for xyz, elements in poses
    ]


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # 1. Redocking: does the top-ranked pose reproduce the crystal pose?
    print(f"{'complex':8} {'top-1 (A)':>10} {'best of 9 (A)':>14}  status")
    all_rmsds = {}
    for pdb_id, name in COMPLEXES.items():
        rmsds = redock_rmsds(pdb_id)
        all_rmsds[pdb_id] = rmsds
        result = evaluate_redocking_success(rmsds, ranked_poses=True)
        print(f"{pdb_id:8} {rmsds[0]:>10.2f} {min(rmsds):>14.2f}  {result['status']}   ({name})")

    # 2. Screening enrichment, together with the redocking of one complex.
    labels, scores, properties, meta = load_docking_csv(SCREENING_CSV)
    report = assess_docking_quality(
        labels=labels,
        scores=scores,
        rmsd_values=all_rmsds["3PTB"],
        rmsd_ranked_poses=True,
        active_properties={k: v[labels == 1] for k, v in properties.items()},
        decoy_properties={k: v[labels == 0] for k, v in properties.items()},
        lower_is_better=True,
    )
    print(f"\nOverall status: {report.overall_status}")
    for item in report.enrichment_metrics.values():
        ci = f"[{item.ci_lower_95:.3f}, {item.ci_upper_95:.3f}]" if item.ci_lower_95 is not None else ""
        print(f" - {item.name:18s} {item.value:8.3f}  {ci:18s} {item.status}")

    # 3. Write the report, figures, tables and methods paragraph.
    generate_docking_figures(labels, scores, report, OUTPUT_DIR, lower_is_better=True,
                             rmsd_values=all_rmsds["3PTB"])
    assets = generate_docking_manuscript_assets(report, OUTPUT_DIR)
    with open(assets["methods_text"], "r", encoding="utf-8") as f:
        methods_text = f.read()
    with open(assets["citation_bib"], "r", encoding="utf-8") as f:
        citation_bib = f.read()
    html_path = os.path.join(OUTPUT_DIR, "report.html")
    generate_docking_html_report(report, html_path, methods_text=methods_text, citation_bib=citation_bib)
    print(f"\nReport written to {html_path}")


if __name__ == "__main__":
    main()
