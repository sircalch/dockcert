"""
Enrichment metrics against the reference implementation in RDKit (rdkit.ML.Scoring), on scores
without ties. With tied scores RDKit's result depends on the input order, whereas DockCert treats
ties symmetrically, so small differences are expected there.
"""
import numpy as np
import pytest
from dockcert.core.enrichment import calculate_roc_auc, calculate_bedroc, calculate_rie, calculate_enrichment_factor

Scoring = pytest.importorskip("rdkit.ML.Scoring.Scoring")


@pytest.mark.parametrize("seed", [0, 1, 2])
def test_metrics_match_rdkit(seed):
    rng = np.random.default_rng(seed)
    n_act, n_dec = 40, 1500
    y = np.r_[np.ones(n_act, int), np.zeros(n_dec, int)]
    s = np.r_[rng.normal(-9.0, 1.0, n_act), rng.normal(-7.0, 1.0, n_dec)]   # lower = better, no ties
    order = np.argsort(s)
    ranked = [[s[i], y[i]] for i in order]
    for a in (20.0, 80.5):
        assert calculate_bedroc(y, s, a) == pytest.approx(Scoring.CalcBEDROC(ranked, 1, a), abs=1e-9)
        assert calculate_rie(y, s, a) == pytest.approx(Scoring.CalcRIE(ranked, 1, a), abs=1e-9)
    assert calculate_roc_auc(y, s) == pytest.approx(Scoring.CalcAUC(ranked, 1), abs=1e-9)
    for f in (0.005, 0.01, 0.05):
        assert calculate_enrichment_factor(y, s, f)[0] == pytest.approx(Scoring.CalcEnrichment(ranked, 1, [f])[0], abs=1e-9)


def test_delong_matches_bruteforce():
    """DeLong variance from ranks equals the direct computation from the full psi matrix (with ties)."""
    from scipy.stats import norm
    from dockcert.core.enrichment import calculate_roc_auc_ci
    rng = np.random.default_rng(3)
    y = np.r_[np.ones(25, int), np.zeros(400, int)]
    s = np.r_[rng.normal(-1.5, 1, 25), np.round(rng.normal(0, 1, 400), 1)]
    a, d = -s[y == 1], -s[y == 0]
    psi = (a[:, None] > d[None, :]) + 0.5 * (a[:, None] == d[None, :])
    auc = psi.mean()
    var = psi.mean(1).var(ddof=1) / len(a) + psi.mean(0).var(ddof=1) / len(d)
    lg, h = np.log(auc / (1 - auc)), norm.ppf(0.975) * np.sqrt(var) / (auc * (1 - auc))
    got = calculate_roc_auc_ci(y, s)
    assert got[0] == pytest.approx(auc)
    assert got[1] == pytest.approx(1 / (1 + np.exp(-(lg - h))))
    assert got[2] == pytest.approx(1 / (1 + np.exp(-(lg + h))))
