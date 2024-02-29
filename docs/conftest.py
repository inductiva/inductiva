import pytest
import inductiva

def pytest_markdown_docs_globals():
    return {"inductiva": inductiva}
