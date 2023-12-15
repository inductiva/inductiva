# pylint: disable=missing-module-docstring
import random
import re
import string


def split_camel_case(s):
    """Split a camel case string into a list of strings."""
    return re.sub("([A-Z][a-z]+)", r" \1", re.sub("([A-Z]+)", r" \1",
                                                  s)).split()


def create_random_tag(size: int = 4):

    samples = [random.choice(string.ascii_lowercase) for _ in range(size)]
    tag = "".join(samples)

    return tag
