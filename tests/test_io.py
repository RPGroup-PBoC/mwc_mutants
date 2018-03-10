import pytest
import numpy as np
import pandas as pd
import sys
sys.path.insert(0, '../')
from mwc.io import scrape_frontmatter


def test_scrape_frontmatter():
    accept_dict = {'status': 'accept', 'reason': 'test reason'}
    reject_dict = {'status': 'reject', 'reason': 'test reason'}
    questionable_dict = {'status': 'questionable',
                         'reason': 'test reason'}
    dicts = {'test_accept.md': accept_dict, 'test_reject.md': reject_dict,
             'test_questionable.md': questionable_dict}

    def cmp_dict(truth, test_file, directory):
        test = scrape_frontmatter(directory, file=test_file)
        assert test.keys() == truth.keys()
        assert test.items() == test.items()

    for k, v in zip(dicts.keys(), dicts.values()):
        cmp_dict(v, k, './test_data', False)
        cmp_dict(v, k, './test_data/', False)

    with pytest.raises(UserWarning):
        assert scrape_frontmatter('./test_data', file='test_wrongstatus.md')
        assert scrape_frontmatter('./test_data', file='test_nostatus.md')
