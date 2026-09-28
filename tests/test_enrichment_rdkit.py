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
