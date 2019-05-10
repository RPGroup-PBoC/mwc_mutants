import numpy as np
import pytest
import sys
sys.path.insert(0, '../')
import mut.thermo

def test_MWC():
    # Bona-fide active probability.
    def pact(c, ka, ki, ep_ai, n):
        numer = (1 + c / ka)**n
        denom = numer + np.exp(-ep_ai) * (1 + c / ki)**n
        return numer / denom

    # Test exceptions.
    with pytest.raises(RuntimeError):
        mwc.thermo.MWC(-1, -1, -1, -1, -1)

    with pytest.raises(ValueError):
        mwc.thermo.MWC(0, 0, 0, 0, 0)

    # Test the pact functions.
    ka, ki, c, ep_ai, n = (100, 1, 100, -1, 2)
    test_pact = pact(c, ka, ki, ep_ai, n)
    assert test_pact == mwc.thermo.MWC(effector_conc=c, ka=ka,
                                       ki=ki, ep_ai=ep_ai, n_sites=n).pact()

    # Test the log transform.
    ep_a, ep_i = 100, -1
    test_pact = pact(c, np.exp(ep_a), np.exp(ep_i), ep_ai, n)
    assert test_pact == mwc.thermo.MWC(
        c, ep_a, ep_i, ep_ai, n, log_transform=True).pact()

    # Test saturation and leakiness.
    leak = (1 + np.exp(-ep_ai))**-1
    sat = (1 + np.exp(-ep_ai) * (ka / ki)**n)**-1
    assert sat == mwc.thermo.MWC(c, ka, ki, ep_ai).saturation()
    assert leak == mwc.thermo.MWC(c, ka, ki, ep_ai).leakiness()

    # Test that everything works with single-arrays.
    c_range = np.logspace(0, 5, 10)
    ka_range = np.logspace(0, 5, 10)
    ki_range = np.logspace(0, 5, 10)
    ep_ai_range = np.logspace(0, 5, 10)
    n_range = np.linspace(0, 10, 10)
    test_c_pact = pact(c_range, ka, ki, ep_ai, n)
    test_ka_pact = pact(c, ka_range, ki, ep_ai, n)
    test_ki_pact = pact(c, ka, ki_range, ep_ai, n)
    test_ep_ai_pact = pact(c, ka, ki, ep_ai_range, n)
    test_n_pact = pact(c, ka, ki, ep_ai, n_range)
    assert (test_c_pact == mwc.thermo.MWC(c_range, ka, ki, ep_ai,
                                          n).pact()).all()
    assert (test_ka_pact == mwc.thermo.MWC(c, ka_range, ki, ep_ai,
                                           n).pact()).all()
    assert (test_ki_pact == mwc.thermo.MWC(c, ka, ki_range, ep_ai,
                                           n).pact()).all()
    assert (test_ep_ai_pact == mwc.thermo.MWC(c, ka, ki, ep_ai_range,
                                              n).pact()).all()
    assert (test_n_pact == mwc.thermo.MWC(c, ka, ki, ep_ai,
                                          n_range).pact()).all()

    # Test meshed variables.
    c, ka, ki, ep_ai, n = np.meshgrid(c_range, ka_range, ki_range,
                                      ep_ai_range, n_range)
    test_pact = pact(c, ka, ki, ep_ai, n)
    assert (test_pact == mwc.thermo.MWC(c, ka, ki, ep_ai, n).pact()).all()

    # Ensure that the values for pact are in the range zero to one.
    max_val, min_val = np.max(test_pact), np.min(test_pact)
    assert (max_val <= 1) & (min_val >= 0)
