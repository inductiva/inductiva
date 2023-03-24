"""Functions that perform dummy operations on the API useful for testing."""
from inductiva.api import invoke_api_from_fn_ptr
# pylint: disable=unused-argument, disable=redefined-builtin


def sleep(secs: float) -> None:
    return invoke_api_from_fn_ptr(locals(), sleep)


# pylint: enable=unused-argument, enable=redefined-builtin
