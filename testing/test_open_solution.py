import pickle
import json
import pytest
import sys
import NGTS_workpackage.wcs_fitting as a
import os
from astropy.io import fits
import numpy as np


@pytest.fixture
def solution(root_dir):
    with open(os.path.join(root_dir, 'testing', 'fixtures', 'wcs_params.json')) as infile:
        return json.load(infile)


@pytest.fixture
def pickle_file(solution, tmpdir):
    fname = str(tmpdir.join('solution.p'))
    with open(fname, 'w') as outfile:
        pickle.dump(solution['wcs'], outfile, protocol=2)
    return fname


@pytest.fixture
def json_file(solution, tmpdir):
    fname = str(tmpdir.join('solution.json'))
    with open(fname, 'w') as outfile:
        json.dump(solution, outfile, indent=2)
    return fname


@pytest.fixture
def fits_file(solution, tmpdir):
    phdu = fits.PrimaryHDU()
    for key in solution['wcs']:
        phdu.header[key] = solution['wcs'][key]

    out_fname = str(tmpdir.join('solution.fits'))
    phdu.writeto(out_fname)
    return out_fname


def assert_same_solutions(a, b):
    nonfloat_keys = ['CTYPE1', 'CTYPE2']
    keys = set(b.keys()) - set(nonfloat_keys)

    for key in nonfloat_keys:
        assert a[key] == b[key]

    assert np.allclose([a[key] for key in keys],
                       [b[key] for key in keys])

def test_load_pickle_file(pickle_file, solution):
    dist_map = a.extract_dist_map(pickle_file)
    assert_same_solutions(dist_map, solution['wcs'])


def test_load_json_file(json_file, solution):
    dist_map = a.extract_dist_map(json_file)
    assert_same_solutions(dist_map, solution['wcs'])


def test_load_fits_file(fits_file, solution):
    dist_map = a.extract_dist_map(fits_file)
    assert_same_solutions(dist_map, solution['wcs'])


def test_exception_for_invalid_filetype():
    with pytest.raises(ValueError) as err:
        a.extract_dist_map('test.bad')
