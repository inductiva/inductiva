"""Configuration for pytest and pytest-markdown."""
import pytest
import inductiva


@pytest.fixture
def simulator():
    fds = inductiva.simulators.FDS()
    return fds


@pytest.fixture
def input_dir():
    input_files_dir = inductiva.utils.download_from_url(
        "https://storage.googleapis.com/inductiva-api-demo-files/"
        "fds-input-example.zip",
        unzip=True)
    return input_files_dir


def pytest_markdown_docs_globals():
    return {"inductiva": inductiva}
