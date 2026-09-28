"""
Enrichment analysis of the DUD-E screens (dude_screen.py) with DockCert, cross-checked against RDKit.

For each target: ROC-AUC with the DeLong (logit) interval, BEDROC (alpha 20 and 80.5), EF at 1% and
5% with stratified bootstrap intervals (1000 resamples), and the same metrics from rdkit.ML.Scoring.
Ligands whose preparation or docking failed are reported and excluded. When a ligand appears more
than once in a score file (a retried run), the last row is used.

Usage: python analyze_dude.py SCORES_DIR OUT_DIR       (SCORES_DIR holds <target>_scores.csv)
"""
import glob
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from dockcert.core.enrichment import (calculate_roc_auc_ci, calculate_bedroc, calculate_enrichment_factor,   # noqa: E402
                                      calculate_roc_auc)
from dockcert.core.bootstrap import bootstrap_enrichment_ci                                                  # noqa: E402


def load(path):
    d = pd.read_csv(path, dtype={"error": str}).drop_duplicates("id", keep="last")
    ok = d[d.vina_score.notna()].copy()
    return ok, d[d.vina_score.isna()]


def rdkit_metrics(y, s):
    from rdkit.ML.Scoring import Scoring
    order = np.argsort(s, kind="stable")
    ranked = [[s[i], y[i]] for i in order]
    return {"auc": Scoring.CalcAUC(ranked, 1), "bedroc20": Scoring.CalcBEDROC(ranked, 1, 20.0),
            "bedroc80": Scoring.CalcBEDROC(ranked, 1, 80.5),
            "ef1": Scoring.CalcEnrichment(ranked, 1, [0.01])[0], "ef5": Scoring.CalcEnrichment(ranked, 1, [0.05])[0]}


def main(scores_dir, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    rows = []
    for f in sorted(glob.glob(os.path.join(scores_dir, "*_scores.csv"))):
        target = os.path.basename(f).split("_scores")[0]
        ok, failed = load(f)
        y = ok.label.astype(int).values
        s = ok.vina_score.astype(float).values
        auc, lo, hi = calculate_roc_auc_ci(y, s)
        ci = bootstrap_enrichment_ci(y, s, n_resamples=1000, random_state=42)
        rd = rdkit_metrics(y, s)
        n_ties = int(len(s) - len(np.unique(s)))
        row = {"target": target, "actives": int(y.sum()), "decoys": int((1 - y).sum()),
               "failed": len(failed), "tied_scores": n_ties,
               "auc": auc, "auc_lo": lo, "auc_hi": hi, "auc_boot_lo": ci["roc_auc"][1], "auc_boot_hi": ci["roc_auc"][2],
               "bedroc20": ci["bedroc_20"][0], "bedroc20_lo": ci["bedroc_20"][1], "bedroc20_hi": ci["bedroc_20"][2],
               "bedroc80": ci["bedroc_80"][0], "bedroc80_lo": ci["bedroc_80"][1], "bedroc80_hi": ci["bedroc_80"][2],
               "ef1": ci["ef_1pct"][0], "ef1_lo": ci["ef_1pct"][1], "ef1_hi": ci["ef_1pct"][2],
               "ef1_max": calculate_enrichment_factor(y, s, 0.01)[1],
               "ef5": ci["ef_5pct"][0], "ef5_lo": ci["ef_5pct"][1], "ef5_hi": ci["ef_5pct"][2]}
        for k, v in rd.items():
            row[f"rdkit_{k}"] = v
        rows.append(row)
        print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in row.items()})
        ok.to_csv(os.path.join(out_dir, f"dude_{target}_scores.csv"), index=False)
    pd.DataFrame(rows).to_csv(os.path.join(out_dir, "dude_metrics.csv"), index=False)


if __name__ == "__main__":
    main(*sys.argv[1:3])
