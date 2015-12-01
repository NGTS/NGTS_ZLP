import pytest
import shutil
import os
import json

from NGTS_workpackage import wcs_fitting as w

TESTDIR = os.path.dirname(__file__)


@pytest.fixture
def source_filename(tmpdir):
    fname = str(tmpdir.join('unsolved-image.fits'))
    shutil.copyfile(os.path.join(TESTDIR, 'fixtures', 'unsolved-image.fits'), fname)
    return fname


@pytest.fixture
def dist_map():
    fname = os.path.join(TESTDIR, 'fixtures', 'wcs_params.json')
    return w.extract_dist_map(fname)


@pytest.fixture
def wcsref():
    return os.path.join(TESTDIR, 'fixtures', 'reference-catalogue.fits')


def test_solve_with_casu_solve(source_filename, dist_map, wcsref):
    w.casu_solve(source_filename, wcsref, dist_map, thresh=7)
    assert os.path.isfile(source_filename.replace('.fits', '') + '.png')
