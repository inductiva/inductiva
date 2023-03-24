"""Functions that perform dummy operations on the API useful for testing."""
from .api import invoke_api
# pylint: disable=unused-argument, disable=redefined-builtin


def sleep(secs: float) -> None:
    return invoke_api(locals(), sleep)


# pylint: enable=unused-argument, enable=redefined-builtin
