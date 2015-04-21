try:
    from unittest import mock
except ImportError:
    import mock

from NGTS_workpackage import wcs_fitting as w


@mock.patch('NGTS_workpackage.catmatch.fits.open')
def test_function_apply_correct(fits_open):
    dicty = {'blah': 10}
    casuin = mock.MagicMock()

    w.apply_correct(dicty, casuin)

    infile = fits_open.return_value.__enter__.return_value
    assert mock.call(0) in infile.__getitem__.call_args_list
    header = infile.__getitem__.return_value.header
    assert mock.call('blah', 10) in header.__setitem__.call_args_list
