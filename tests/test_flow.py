import fcsparser
import pandas as pd
import numpy as np
import pytest
import scipy.stats
import sys
sys.path.insert(0, '../')
from mwc.flow import *
np.random.seed(42)
test_path = fcsparser.test_sample_path
metadata, data = fcsparser.parse(test_path)


def test_gaussian_gate():

    # Define a MV gaussian pdf.
    mean = [1E4, 1E4]
    cov = [[100, 0], [0, 100]]
    rv = scipy.stats.multivariate_normal(mean=mean, cov=cov)

    # Draw random samples
    n_samples = 100000
    rnd = rv.rvs(n_samples)
    rnd_df = pd.DataFrame(rnd, columns=['x', 'y'])
    rnd_df['pdf'] = rv.pdf(rnd_df)
    rnd_df.sort_values(by='pdf', inplace=True)

    # Identify the percentile.
    alpha_range = [0.01, 0.1, 0.2, 0.4, 0.6, 0.8, 1]

    # Gate for a range of alpha
    def _test_fn(alpha_range, log_bool):
        for alpha in alpha_range:
            fn_gated = gaussian_gate(rnd_df, alpha,  x_val='x', y_val='y',
                                     log=log_bool)
            frac = len(fn_gated) / n_samples
            assert alpha == pytest.approx(frac, 0.1)

            # Gate the fake data by percentile.
            bound = np.percentile(rnd_df['pdf'], (1 - alpha) * 100)
            gated = rnd_df.loc[rnd_df['pdf'] > bound]

            # Compute mean and std for means
            mean_perc = gated[['x', 'y']].mean().values
            std_perc = gated[['x', 'y']].std().values
            mean_fn = fn_gated[['x', 'y']].mean().values
            std_fn = fn_gated[['x', 'y']].std().values

            # Ensure the computed means and std dev. are within 0.15 units
            assert mean_perc == pytest.approx(mean_fn, 0.1)
            assert std_perc == pytest.approx(std_fn, 0.1)

    # Run the tests
    _test_fn(alpha_range, False)
    _test_fn(alpha_range, True)
