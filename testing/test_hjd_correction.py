from NGTS_workpackage.hjd_correction import (
    compute_hjd, compute_hjd_column, append_hjd_column
)
import os
import pytest
import fitsio
import shutil
from numpy import allclose, dtype, float64

ref_jd = 2456834.3842939814
ref_ra = 210.150833
ref_dec = -60.002286
sun_ra = 94.6073456
sun_dec = 23.365036


def get_column_names(fname):
    with fitsio.FITS(fname) as infile:
        hdu = infile[1]
        return hdu.get_colnames()


@pytest.fixture
def example_file():
    return os.path.join(
        os.path.dirname(__file__),
        'hjd_example_catalogue.fits'
    )


@pytest.fixture
def backup_file(example_file, tmpdir):
    out_dir = str(tmpdir)
    out_path = os.path.join(out_dir, os.path.basename(example_file))
    shutil.copyfile(example_file, out_path)
    return out_path


def test_hjd_computation():
    result = compute_hjd(ref_jd, ref_ra, ref_dec, sun_ra, sun_dec)
    expected = 2456834.38785
    assert allclose(result, expected)


def test_read_from_file(example_file):
    hjd_column_data = compute_hjd_column(example_file)
    assert allclose(hjd_column_data[0], 56833.8861)


def test_update_file_column_present(backup_file):
    before_column_names = get_column_names(backup_file)
    assert 'hjd' not in before_column_names

    append_hjd_column(backup_file, column_name='hjd')

    after_column_names = get_column_names(backup_file)
    assert 'hjd' in after_column_names


def test_update_file_column_dtype(backup_file):
    append_hjd_column(backup_file, column_name='hjd')

    with fitsio.FITS(backup_file) as infile:
        column = infile[1]['hjd'].read()

    assert column.dtype.type == dtype(float64)


def test_update_file_values_correct(backup_file):
    append_hjd_column(backup_file)

    with fitsio.FITS(backup_file) as infile:
        catalogue = infile[1]
        ra, dec = [catalogue[key].read() for key in ['ra', 'dec']]
        hjd_column_data = catalogue['hjd'].read()

        header = catalogue.read_header()

    sun_ra, sun_dec = [header[key] for key in ['sun_ra', 'sun_dec']]
    mjd = header['mjd']

    hjd_values = compute_hjd(mjd, ra, dec, sun_ra, sun_dec)
    assert allclose(hjd_column_data, hjd_values, atol=0, rtol=1E-12)
