"""
Known-answer test of DockCert's stratified bootstrap confidence intervals for enrichment metrics.

Binormal model (lower score = better): decoy scores ~ N(0, 1), active scores ~ N(-mu, 1). For this
model the population values are known:
  ROC-AUC   = Phi(mu / sqrt(2))
  EF_x      = F_a(t) / x, where t solves r F_a(t) + (1 - r) F_d(t) = x and r is the active fraction
            (the limit of the enrichment factor for an infinitely large library with the same r)
  BEDROC_a  population value estimated once from a library 1000 times larger (same r)
For each design (actives A, decoys 50 A, separation mu) we draw R independent libraries, compute
DockCert's 95% bootstrap interval (1000 stratified resamples) and count how often it contains the
population value. The Monte Carlo 95% range of the coverage for R replicates is 0.95 +/- 1.96
sqrt(0.95 0.05 / R).

Usage: python bootstrap_coverage.py OUT_CSV [--reps 300] [--workers 8]
"""
import argparse
import os
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd
from scipy.optimize import brentq
from scipy.stats import norm

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from dockcert.core.bootstrap import bootstrap_enrichment_ci      # noqa: E402
from dockcert.core.enrichment import calculate_bedroc            # noqa: E402

KEYS = {"roc_auc": "auc", "ef_1pct": "ef1", "ef_5pct": "ef5", "bedroc_20": "bedroc20"}


def population_values(mu, ratio, seed=7):
    r = ratio
    truth = {"auc": float(norm.cdf(mu / np.sqrt(2.0)))}
    for x, key in ((0.01, "ef1"), (0.05, "ef5")):
        t = brentq(lambda t: r * norm.cdf(t + mu) + (1 - r) * norm.cdf(t) - x, -20, 20)
        truth[key] = float(norm.cdf(t + mu) / x)
    rng = np.random.default_rng(seed)
    n_a = 50_000
    n_d = int(round(n_a * (1 - r) / r))
    y = np.r_[np.ones(n_a, int), np.zeros(n_d, int)]
    s = np.r_[rng.normal(-mu, 1, n_a), rng.normal(0, 1, n_d)]
    truth["bedroc20"] = calculate_bedroc(y, s, alpha=20.0)
    return truth


def one_rep(job):
    n_act, n_dec, mu, seed = job
    rng = np.random.default_rng(seed)
    y = np.r_[np.ones(n_act, int), np.zeros(n_dec, int)]
    s = np.r_[rng.normal(-mu, 1, n_act), rng.normal(0, 1, n_dec)]
    ci = bootstrap_enrichment_ci(y, s, n_resamples=1000, random_state=seed)
    return {v: ci[k] for k, v in KEYS.items() if k in ci}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out_csv")
    ap.add_argument("--reps", type=int, default=300)
    ap.add_argument("--workers", type=int, default=8)
    a = ap.parse_args()
    rows = []
    for mu in (1.0, 2.0):
        for n_act in (20, 50, 100):
            n_dec = 50 * n_act
            truth = population_values(mu, n_act / (n_act + n_dec))
            jobs = [(n_act, n_dec, mu, 1000 * n_act + k) for k in range(a.reps)]
            with ProcessPoolExecutor(max_workers=a.workers) as ex:
                res = list(ex.map(one_rep, jobs, chunksize=4))
            for key, tv in truth.items():
                est = np.array([r[key] for r in res])            # (point, low, high)
                cover = np.mean((est[:, 1] <= tv) & (tv <= est[:, 2]))
                rows.append({"mu": mu, "actives": n_act, "decoys": n_dec, "metric": key, "true": tv,
                             "mean_estimate": est[:, 0].mean(), "bias": est[:, 0].mean() - tv,
                             "coverage": cover, "median_width": float(np.median(est[:, 2] - est[:, 1])),
                             "reps": a.reps})
            print(pd.DataFrame(rows).tail(len(truth)).to_string(index=False), flush=True)
    pd.DataFrame(rows).to_csv(a.out_csv, index=False)


if __name__ == "__main__":
    main()
