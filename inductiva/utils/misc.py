# pylint: disable=missing-module-docstring
import glob
import re


def split_camel_case(s):
    """Split a camel case string into a list of strings."""
    return re.sub("([A-Z][a-z]+)", r" \1", re.sub("([A-Z]+)", r" \1",
                                                  s)).split()


def gen_suffix(name, filenames):
    """Generate a suffix for a filename, based on a list of files."""
    regex = f"^{name}#([0-9]+$)"

    max_suffix = -1
    for filename in filenames:
        if match := re.match(regex, filename):
            max_suffix = max(max_suffix, int(match.group(1)))

    exists_solo = name in filenames
    exists_derived = max_suffix >= 0

    if exists_derived:
        suffix = f"#{max_suffix+1}"
    elif exists_solo:
        suffix = "#2"
    else:
        suffix = ""
    return suffix


def gen_name(name):
    name = name.strip()
    filenames = glob.glob(f"{name}*")
    return name + gen_suffix(name, filenames)
