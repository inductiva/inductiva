# pylint: disable=missing-module-docstring
import re


def split_camel_case(s):
    """Split a camel case string into a list of strings."""
    return re.sub("([A-Z][a-z]+)", r" \1", re.sub("([A-Z]+)", r" \1",
                                                  s)).split()
