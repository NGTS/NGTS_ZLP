import pytest


@pytest.fixture
def source_filename(tmpdir):
    return 'test'


@pytest.fixture
def dist_map():
    return 'dist_map.p'


def test_solve_with_casu_solve(source_filename, dist_map):
    pass
