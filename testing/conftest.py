import pytest
import os


@pytest.fixture
def root_dir():
    return os.path.realpath(os.path.join(os.path.dirname(__file__), '..'))
