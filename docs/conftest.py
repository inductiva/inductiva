"""Configuration for pytest and pytest-markdown."""
import pytest
import inductiva


@pytest.fixture
def hello_world():
    return "Hello, World!"


def pytest_markdown_docs_globals():
    return {"inductiva": inductiva}
