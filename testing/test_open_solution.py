import pickle
import json
import pytest
import sys
import NGTS_workpackage.wcs_fitting as a
import os


@pytest.fixture
def solution(root_dir):
    with open(os.path.join(root_dir, 'testing', 'fixtures', 'wcs_params.json')) as infile:
        return json.load(infile)['wcs']


@pytest.fixture
def pickle_file(solution, tmpdir):
    fname = str(tmpdir.join('solution.p'))
    with open(fname, 'w') as outfile:
        pickle.dump(solution, outfile, protocol=2)
    return fname


@pytest.fixture
def json_file(solution, tmpdir):
    fname = str(tmpdir.join('solution.p'))
    with open(fname, 'w') as outfile:
        json.dump(solution, outfile, indent=2)
    return fname


def test_load_pickle_file(pickle_file, solution):
    dist_map = a.extract_dist_map(pickle_file)
    assert dist_map == solution


def test_load_json_file(json_file, solution):
    dist_map = a.extract_dist_map(json_file)
    assert dist_map == solution
