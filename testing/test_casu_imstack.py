import mock
from NGTS_workpackage import casutools as casu

@mock.patch('NGTS_workpackage.casutools.run_command')
def test_imstack(mock_run_command):
    print mock_run_command
    assert False

