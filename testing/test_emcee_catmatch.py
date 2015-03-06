import pytest
import mock
import sys
sys.path.insert(0, 'bin')
import emcee_catmatch


def test_extract_from_index_file():
    lines = ['1 2 -2.0 2.0 -4.0 4.0', ]
    expected = (
            ['1', ],
            [[-2., 2.]],
            [[-4., 4.]],
            )

    with mock.patch('emcee_catmatch.open', create=True) as m:
        m.return_value = mock.MagicMock(spec=file)
        handle = m.return_value.__enter__.return_value
        handle.__iter__.return_value = lines
        assert emcee_catmatch.extract_coordinate_limits(None) == expected
