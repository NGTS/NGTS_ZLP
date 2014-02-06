import mock
from NGTS_workpackage import casutools
import os

class printer(object):
    def __init__(self):
        self.args = []

    def __call__(self, *args, **kwargs):
        self.args.append(["Args: {}".format(args), "kwargs: {}".format(kwargs)])

@mock.patch('NGTS_workpackage.casutools.run_command', new_callable=printer)
def test_imstack(mock_run_command):
    filelist = 'files.txt'
    confidence_map = 'confidence.fits'

    casutools.imstack(filelist, confidence_map)
    assert len(mock_run_command.args) == 1

    command = mock_run_command.args[0]
    assert filelist in command[0] and confidence_map in command[0]

def test_construct_filelist():
    filelist = './files.txt'
    assert casutools.construct_filelist_argument(filelist) == '@{}'.format(filelist)

def run_command():
    pwd = os.getcwd()
    assert casutools.run_command(['pwd']) == pwd
