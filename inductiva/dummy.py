"""
Functions that call `dummy` endpoints.
"""
from .api import invoke_api

# pylint: disable=unused-argument, disable=redefined-builtin


def sum(a: float, b: float) -> float:
    return invoke_api(locals(), sum)


# TODO: this one creates a file. How to specify it as output?
#def gen_file(size: int, sleep_s: int) -> float:
#   return invoke_api(locals(), gen_file)

# pylint: enable=unused-argument, enable=redefined-builtin
