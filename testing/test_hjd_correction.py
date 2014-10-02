from NGTS_workpackage.hjd_correction import (
    compute_hjd_correction, compute_hjd_correction_column,
    append_hjd_correction_column
)
import os
import pytest
import fitsio
import shutil
from numpy import allclose, dtype, float64


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
    ref_jd = 2456834.4
    ref_ra = 14. * 15.
    ref_dec = -60.0
    sun_ra = 6.30860 * 15.
    sun_dec = 23.37027

    result = compute_hjd_correction(ref_jd, ref_ra, ref_dec, sun_ra, sun_dec)
    expected = 2456834.403118045 - 2456834.399988426
    assert allclose(result, expected, rtol=1E-2, atol=0)


def test_read_from_file(example_file):
    hjd_column_data = compute_hjd_correction_column(example_file)
    assert allclose(hjd_column_data[0], 0.0018173871387605742)


def test_update_file_column_present(backup_file):
    before_column_names = get_column_names(backup_file)
    assert 'hjd' not in before_column_names

    append_hjd_correction_column(backup_file, column_name='hjd')

    after_column_names = get_column_names(backup_file)
    assert 'hjd' in after_column_names


def test_default_column_name(backup_file):
    target = 'hjd_correction'
    before_column_names = get_column_names(backup_file)
    assert target not in before_column_names

    append_hjd_correction_column(backup_file)

    after_column_names = get_column_names(backup_file)
    assert target in after_column_names and 'hjd' not in after_column_names


def test_update_file_column_dtype(backup_file):
    append_hjd_correction_column(backup_file, column_name='hjd_correction')

    with fitsio.FITS(backup_file) as infile:
        column = infile[1]['hjd_correction'].read()

    assert column.dtype.type == dtype(float64)


def test_update_file_values_correct(backup_file):
    append_hjd_correction_column(backup_file, column_name='hjd_correction')

    with fitsio.FITS(backup_file) as infile:
        catalogue = infile[1]
        ra, dec = [catalogue[key].read() for key in ['ra', 'dec']]
        hjd_column_data = catalogue['hjd_correction'].read()

        header = catalogue.read_header()

    sun_ra, sun_dec = [header[key] for key in ['sun_ra', 'sun_dec']]
    mjd = header['mjd']

    hjd_values = compute_hjd_correction(mjd, ra, dec, sun_ra, sun_dec)
    assert allclose(hjd_column_data, hjd_values)
