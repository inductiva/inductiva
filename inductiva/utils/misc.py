# pylint: disable=missing-module-docstring
import re


def split_camel_case(s):
    """Split a camel case string into a list of strings."""
    return re.sub("([A-Z][a-z]+)", r" \1", re.sub("([A-Z]+)", r" \1",
                                                  s)).split()


def format_bytes(n_bytes: int) -> str:
    """Convert bytes to human readable string."""
    res = float(n_bytes)

    for unit in ["bytes", "KiB", "MiB", "GiB", "TiB"]:
        if res < 1024:
            return f"{res:.2f} {unit}"
        res /= 1024

    return f"{format_bytes:.2f} PiB"
