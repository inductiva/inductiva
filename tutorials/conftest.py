"""Configuration for pytest and pytest-markdown."""
import pytest
import datetime
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


@pytest.fixture(scope="session")
def machine_group_sample():
    mg = inductiva.resources.MachineGroup(
        "e2-micro",
        max_idle_time=datetime.timedelta(minutes=1),
    )
    mg.start()
    yield mg
    mg.terminate()
