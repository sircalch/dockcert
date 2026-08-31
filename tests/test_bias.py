"""
Tests for decoy bias evaluation and physicochemical property distribution checks.
"""

import numpy as np
import pytest
from dockcert.core.bias import evaluate_decoy_bias


def test_decoy_bias_unbiased():
    rng = np.random.default_rng(42)
    # Identical MW distributions
    act_mw = rng.normal(400, 30, size=100)
    dec_mw = rng.normal(400, 30, size=1000)
    
    res = evaluate_decoy_bias({"MW": act_mw}, {"MW": dec_mw})
    assert res["status"] == "PASS"
    assert res["risk_level"] == "LOW"


def test_decoy_bias_severe():
    rng = np.random.default_rng(42)
    # Actives are 500 Da, decoys are 250 Da
    act_mw = rng.normal(500, 20, size=100)
    dec_mw = rng.normal(250, 20, size=1000)
    
    act_logp = rng.normal(4.5, 0.3, size=100)
    dec_logp = rng.normal(1.2, 0.3, size=1000)
    
    res = evaluate_decoy_bias(
        {"MW": act_mw, "LogP": act_logp, "HBD": act_mw*0.01},
        {"MW": dec_mw, "LogP": dec_logp, "HBD": dec_mw*0.01}
    )
    assert res["risk_level"] == "HIGH"
    assert res["status"] == "WARNING"
