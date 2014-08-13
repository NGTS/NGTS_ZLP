import sys
sys.path.extend(['.', 'testdata'])
import fitsio
import pytest
import numpy as np

from NGTS_workpackage.wcs_status import wcs_succeeded, set_wcs_status

@pytest.fixture
def blank_fits(tmpdir):
    fname = str(tmpdir.join("test.fits"))
    with fitsio.FITS(fname, 'rw', clobber=True) as outfile:
        outfile.write(np.zeros((2, 2)))

    return fname


@pytest.fixture
def failed_fits_file(blank_fits):
    with fitsio.FITS(blank_fits, 'rw') as outfile:
        outfile[0].write_key('wcscompl', False, comment='WCS succeeded?')

    return blank_fits

@pytest.fixture
def succeeded_fits_file(blank_fits):
    with fitsio.FITS(blank_fits, 'rw') as outfile:
        outfile[0].write_key('wcscompl', True, comment='WCS succeeded?')

    return blank_fits


def test_succeeded_fits_file(succeeded_fits_file):
    assert wcs_succeeded(succeeded_fits_file)

def test_failed_fits_file(failed_fits_file):
    assert not wcs_succeeded(failed_fits_file)

def test_file_without_key(blank_fits):
    header = fitsio.read_header(blank_fits)
    assert 'wcscompl' not in header
    assert wcs_succeeded(blank_fits)

def test_setting_key_failed(blank_fits):
    old_header = fitsio.read_header(blank_fits)
    assert 'wcscompl' not in old_header
    set_wcs_status(blank_fits, succeeded=False)
    assert not wcs_succeeded(blank_fits)

def test_setting_key_succeeded(blank_fits):
    old_header = fitsio.read_header(blank_fits)
    assert 'wcscompl' not in old_header
    set_wcs_status(blank_fits, succeeded=True)
    assert wcs_succeeded(blank_fits)

def test_second_hdu(blank_fits):
    with fitsio.FITS(blank_fits, 'rw') as outfile:
        outfile.write(np.zeros((2, 2)))

        assert len(outfile) == 2

    set_wcs_status(blank_fits, succeeded=True, hdu=1)

    with fitsio.FITS(blank_fits) as infile:
        hdu = infile[1]
        header = hdu.read_header()

        assert header['wcscompl'] == True

    assert wcs_succeeded(blank_fits, hdu=1)
