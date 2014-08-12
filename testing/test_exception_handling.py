import mock
import pytest
import sys
sys.path.extend(['.', 'testdata'])
from NGTS_workpackage.wcs_fitting import handle_errors_in_casu_solve

@mock.patch('NGTS_workpackage.wcs_fitting.casu_solve')
def test_exception_handling(mock_casu_solve):
    casuin = mock.MagicMock()
    wcsref = mock.MagicMock()
    mock_casu_solve.side_effect = RuntimeError("Test exception")

    with pytest.raises(RuntimeError) as err:
        mock_casu_solve()

    handle_errors_in_casu_solve(casuin, wcsref)

