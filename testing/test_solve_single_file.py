try:
    from unittest import mock
except ImportError:
    import mock

from NGTS_workpackage import wcs_fitting as w


@mock.patch('NGTS_workpackage.wcs_fitting.fitsio.FITS')
def test_function_apply_correct(FITS):
    dicty = {'blah': 10}
    casuin = mock.MagicMock()

    w.apply_correct(dicty, casuin)

    fn = FITS.return_value.__enter__.return_value.__getitem__.return_value.write_key
    assert mock.call('blah', 10) in fn.call_args_list


@mock.patch('NGTS_workpackage.wcs_fitting.fitsio.FITS')
def test_apply_correct_does_not_stamp_shifts(FITS):
    dicty = {}
    casuin = mock.MagicMock()

    w.apply_correct(dicty, casuin)

    fn = FITS.return_value.__enter__.return_value.__getitem__.return_value.write_key
    assert fn.call_args_list == []
