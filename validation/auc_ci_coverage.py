"""
Coverage of the DeLong (logit) ROC-AUC interval on the same binormal libraries (same seeds) as
bootstrap_coverage.py, so that the two intervals are compared on identical data.

Usage: python auc_ci_coverage.py OUT_CSV [--reps 300]
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import norm

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from dockcert.core.enrichment import calculate_roc_auc_ci      # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out_csv")
    ap.add_argument("--reps", type=int, default=300)
    a = ap.parse_args()
    rows = []
    for mu in (1.0, 2.0):
        truth = float(norm.cdf(mu / np.sqrt(2.0)))
        for n_act in (20, 50, 100):
            n_dec = 50 * n_act
            hits, widths = [], []
            for k in range(a.reps):
                rng = np.random.default_rng(1000 * n_act + k)          # same draws as bootstrap_coverage
                y = np.r_[np.ones(n_act, int), np.zeros(n_dec, int)]
                s = np.r_[rng.normal(-mu, 1, n_act), rng.normal(0, 1, n_dec)]
                auc, lo, hi = calculate_roc_auc_ci(y, s)
                hits.append(lo <= truth <= hi)
                widths.append(hi - lo)
            rows.append({"mu": mu, "actives": n_act, "decoys": n_dec, "metric": "auc_delong_logit",
                         "true": truth, "coverage": float(np.mean(hits)),
                         "median_width": float(np.median(widths)), "reps": a.reps})
            print(rows[-1], flush=True)
    pd.DataFrame(rows).to_csv(a.out_csv, index=False)


if __name__ == "__main__":
    main()
