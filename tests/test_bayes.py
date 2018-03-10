import numpy as np
import pandas as pd
import scipy.optimize
import statsmodels.tools.numdiff as smnd
import pytest
import sys
sys.path.insert(0, '../')
from mwc.bayes import deterministic_log_posterior, log_prior, estimate_calibration_factor

# Seed the rng for reproducibility
np.random.seed(666)

# Set up a dilution experiment.
alpha = 100
num_sim = 10
n_tot = np.arange(10, 1E3, 1).astype(int)
n1 = np.random.binomial(n_tot, p=0.5)
n2 = n_tot - n1
i1 = alpha * n1
i2 = alpha * n2


def test_estimate_calibration_factor():
    alpha_opt, alpha_std = estimate_calibration_factor(i1, i2)
    assert alpha_opt == pytest.approx(alpha, 0.5)
    assert alpha_std < 0.05 * alpha

    with pytest.raises(ValueError):
        assert estimate_calibration_factor(-i1, i2)
        assert estimate_calibration_factor(i1, -i2)
        assert estimate_calibration_factor(i1, i2, p=1.01)
        assert estimate_calibration_factor(i1, i2, p=-0.01)
