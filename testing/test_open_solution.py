import pickle
import json
import pytest
import sys
sys.path.insert(0, 'bin')
import assess_astrometric_solution as a


@pytest.fixture
def solution():
    return {'a': {'b': 10}}


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
