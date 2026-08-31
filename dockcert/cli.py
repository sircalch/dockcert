"""
Command Line Interface (CLI) for DockCert.
"""

import sys
import os
import argparse
import numpy as np

from dockcert import __version__
from dockcert.parsers.generic_csv import load_docking_csv
from dockcert.parsers.structure_io import load_molecule_coordinates
from dockcert.core.rmsd import calculate_heavy_atom_rmsd, calculate_symmetry_corrected_rmsd
from dockcert.core.scoring import assess_docking_quality
from dockcert.reporters.plot_generator import generate_docking_figures
from dockcert.reporters.manuscript_prep import generate_docking_manuscript_assets
from dockcert.reporters.html_report import generate_docking_html_report


def print_banner():
    banner = rf"""
  _____             _    _____          _   
 |  __ \           | |  / ____|        | |  
 | |  | | ___   ___| | _| |     ___ _ __| |_ 
 | |  | |/ _ \ / __| |/ / |    / _ \ '__| __|
 | |__| | (_) | (__|   <| |___|  __/ |  | |_ 
 |_____/ \___/ \___|_|\_\\_____\___|_|   \__| v{__version__}

 Molecular Docking Validation & Statistical Quality Toolkit
 Monreal-Hernández et al., 2026
"""
    print(banner)


def run_demo(output_dir: str = "dockcert_demo_output"):
    """
    Generates a realistic DUD-E-like benchmark dataset (100 actives, 2000 property-matched decoys,
    redocking RMSD = 1.43 A) and executes the full validation pipeline.
    """
    print(f"\n[DockCert] Running demonstration mode...")
    os.makedirs(output_dir, exist_ok=True)
    
    n_actives = 100
    n_decoys = 2000
    
    rng = np.random.default_rng(42)
    
    # Active scores: N(-9.4 kcal/mol, 1.1)
    active_scores = rng.normal(-9.4, 1.1, size=n_actives)
    # Decoy scores: N(-6.8 kcal/mol, 1.2)
    decoy_scores = rng.normal(-6.8, 1.2, size=n_decoys)
    
    labels = np.concatenate([np.ones(n_actives, dtype=int), np.zeros(n_decoys, dtype=int)])
    scores = np.concatenate([active_scores, decoy_scores])
    
    # Redocking poses: Best RMSD = 1.43 A, ensemble of 9 poses
    redocking_rmsds = [1.43, 1.78, 2.15, 2.40, 2.89, 3.10, 3.45, 4.12, 4.55]
    
    # Physicochemical properties for bias audit (MW & LogP)
    active_mw = rng.normal(380.0, 45.0, size=n_actives)
    decoy_mw = rng.normal(375.0, 50.0, size=n_decoys)
    
    active_logp = rng.normal(2.8, 0.7, size=n_actives)
    decoy_logp = rng.normal(2.7, 0.8, size=n_decoys)
    
    active_props = {"MolecularWeight": active_mw, "LogP": active_logp}
    decoy_props = {"MolecularWeight": decoy_mw, "LogP": decoy_logp}
    
    print("  -> Calculating ROC-AUC, PR-AUC, BEDROC (alpha=20.0), EF1%, EF5%, logAUC, and Bootstrap 95% CIs...")
    report = assess_docking_quality(
        labels=labels,
        scores=scores,
        rmsd_values=redocking_rmsds,
        active_properties=active_props,
        decoy_properties=decoy_props,
        lower_is_better=True
    )
    
    print("  -> Generating publication-quality vector charts (ROC, PR, Score distributions, RMSD)...")
    generate_docking_figures(labels, scores, report, output_dir, lower_is_better=True, rmsd_values=redocking_rmsds)
    
    print("  -> Drafting manuscript Methods text snippet, summary LaTeX tables, and BibTeX citations...")
    assets = generate_docking_manuscript_assets(report, output_dir)
    
    with open(assets["methods_text"], "r", encoding="utf-8") as f:
        methods_txt = f.read()
    with open(assets["citation_bib"], "r", encoding="utf-8") as f:
        bib_txt = f.read()
        
    html_p = os.path.join(output_dir, "report.html")
    print(f"  -> Writing interactive dashboard to {html_p}...")
    generate_docking_html_report(report, html_p, methods_text=methods_txt, citation_bib=bib_txt)
    
    print("\n" + "="*70)
    print(f" [RESULT] Overall Docking Certification: {report.overall_status}")
    print(f" [SCORE]  {report.validation_score}")
    print("="*70)
    if report.redocking_result:
        rr = report.redocking_result
        print(f" * Redocking RMSD     : Best = {rr['min_rmsd']:.2f} A | Success Rate (<=2A) = {rr['success_rate_2a']:.1f}% | Status: {rr['status']}")
    for k, item in report.enrichment_metrics.items():
        ci_s = f"[{item.ci_lower_95:.2f}, {item.ci_upper_95:.2f}]" if item.ci_lower_95 is not None else ""
        print(f" * {item.name:18s}: Value = {item.value:6.3f} {ci_s:14s} | Threshold = {item.threshold_pass:4.2f} | Status: {item.status}")
    if report.decoy_bias_result:
        print(f" * Decoy Bias Risk    : Level = {report.decoy_bias_result['risk_level']} | Status: {report.decoy_bias_result['status']}")
    print("="*70)
    print(f"\nAll outputs successfully saved to: {os.path.abspath(output_dir)}/")
    print(f"Open {os.path.abspath(html_p)} in your browser to inspect the full report.\n")


def run_assess(args):
    """
    Evaluates user-provided CSV and optional structure files.
    """
    output_dir = args.output
    os.makedirs(output_dir, exist_ok=True)
    
    csv_file = args.input
    if not csv_file:
        print("[Error] Please specify a results CSV file with --input.", file=sys.stderr)
        sys.exit(1)
        
    print(f"\n[DockCert] Loading docking screening dataset from {csv_file}...")
    labels, scores, props_dict, meta = load_docking_csv(
        csv_file,
        score_column=args.score_col,
        label_column=args.label_col
    )
    print(f"  -> Identified {meta['n_total']} entries: {meta['n_actives']} actives, {meta['n_decoys']} decoys.")
    
    # RMSD from structures if provided
    rmsd_values = None
    if args.ref_ligand and args.docked_pose:
        print("  -> Calculating heavy-atom RMSD between reference and docked pose...")
        c_ref, el_ref = load_molecule_coordinates(args.ref_ligand)
        c_dock, el_dock = load_molecule_coordinates(args.docked_pose)
        rmsd_val = calculate_symmetry_corrected_rmsd(c_ref, c_dock, elements=el_ref)
        rmsd_values = [rmsd_val]
        print(f"  -> Calculated Redocking RMSD: {rmsd_val:.2f} A")
        
    print("  -> Performing statistical validation and enrichment analysis...")
    report = assess_docking_quality(
        labels=labels,
        scores=scores,
        rmsd_values=rmsd_values,
        lower_is_better=not args.higher_is_better
    )
    
    print("  -> Generating publication charts...")
    generate_docking_figures(labels, scores, report, output_dir, lower_is_better=not args.higher_is_better, rmsd_values=rmsd_values)
    
    print("  -> Generating manuscript text, LaTeX summary table, and BibTeX citations...")
    assets = generate_docking_manuscript_assets(report, output_dir)
    
    with open(assets["methods_text"], "r", encoding="utf-8") as f:
        methods_txt = f.read()
    with open(assets["citation_bib"], "r", encoding="utf-8") as f:
        bib_txt = f.read()
        
    html_p = os.path.join(output_dir, "report.html")
    print(f"  -> Writing HTML quality report to {html_p}...")
    generate_docking_html_report(report, html_p, methods_text=methods_txt, citation_bib=bib_txt)
    
    print("\n" + "="*70)
    print(f" [RESULT] Overall Docking Certification: {report.overall_status}")
    print(f" [SCORE]  {report.validation_score}")
    print("="*70)
    for k, item in report.enrichment_metrics.items():
        print(f" * {item.name:18s}: {item.value:6.3f} | Status: {item.status}")
    if report.redocking_result:
        print(f" * Redocking RMSD     : {report.redocking_result['min_rmsd']:.2f} A | Status: {report.redocking_result['status']}")
    print("="*70)
    print(f"\nReport ready at: {os.path.abspath(html_p)}\n")


def print_citation():
    bib = """@software{monreal2026dockcert,
  author = {Monreal-Hern\\'andez, Andre},
  title = {{DockCert: An Open-Source Toolkit for Statistical Validation, Enrichment Metrics, and Reproducibility Assessment of Molecular Docking Studies}},
  year = {2026},
  version = {1.0.0},
  publisher = {Zenodo},
  url = {https://github.com/amonreal/dockcert}
}"""
    print("\nIf you use DockCert in your publications, please cite:\n")
    print("APA Style:")
    print("Monreal-Hernández, A. (2026). DockCert: An Open-Source Toolkit for Statistical Validation, Enrichment Metrics, and Reproducibility Assessment of Molecular Docking Studies (v1.0.0). Zenodo. https://github.com/amonreal/dockcert\n")
    print("BibTeX:")
    print(bib)
    print()


def main():
    parser = argparse.ArgumentParser(
        prog="dockcert",
        description="DockCert: Automated Statistical Validation, Enrichment Metrics, and Reproducibility Toolkit for Molecular Docking."
    )
    parser.add_argument("-v", "--version", action="version", version=f"dockcert {__version__}")
    
    subparsers = parser.add_subparsers(dest="command", help="Available subcommands")
    
    # Assess command
    assess_parser = subparsers.add_parser("assess", help="Assess virtual screening enrichment and docking accuracy")
    assess_parser.add_argument("-i", "--input", required=True, help="Virtual screening results file (.csv, .tsv, .txt)")
    assess_parser.add_argument("--score-col", help="Column name containing docking scores/affinities")
    assess_parser.add_argument("--label-col", help="Column name containing active/decoy labels")
    assess_parser.add_argument("--ref-ligand", help="Reference crystallographic ligand structure (.sdf, .pdb, .pdbqt)")
    assess_parser.add_argument("--docked-pose", help="Docked pose structure (.sdf, .pdb, .pdbqt)")
    assess_parser.add_argument("--higher-is-better", action="store_true", help="Set if higher score values indicate better affinity")
    assess_parser.add_argument("-o", "--output", default="dockcert_output", help="Directory for output report and assets (default: dockcert_output)")
    
    # Demo command
    demo_parser = subparsers.add_parser("demo", help="Run DockCert on a benchmark DUD-E-like simulation dataset")
    demo_parser.add_argument("-o", "--output", default="dockcert_demo_output", help="Output directory (default: dockcert_demo_output)")
    
    # Cite command
    subparsers.add_parser("cite", help="Display BibTeX and APA citation details")
    
    if len(sys.argv) == 1:
        print_banner()
        parser.print_help()
        sys.exit(0)
        
    args = parser.parse_args()
    
    if args.command == "assess":
        print_banner()
        run_assess(args)
    elif args.command == "demo":
        print_banner()
        run_demo(args.output)
    elif args.command == "cite":
        print_banner()
        print_citation()
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
