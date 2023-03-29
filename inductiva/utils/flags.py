"""Utilities for parsing script flags."""

from typing import List, Union


def cast_list_to_int(
        list_of_strings: Union[List[str], None]) -> Union[List[int], None]:
    """Converts a list of strings into a list of integers.

    Args:
        list_of_strings: List containing string variables.

    Returns:
        List of integers or None, if the list of strings is None.
    """
    if list_of_strings is None:
        return list_of_strings

    list_of_ints = [int(string) for string in list_of_strings]
    return list_of_ints


def cast_list_to_float(
        list_of_strings: Union[List[str], None]) -> Union[List[float], None]:
    """Converts a list of strings into a list of floats.

    Args:
        list_of_strings: List containing string variables.

    Returns:
        List of floats or None, if the list of strings is None.
    """
    if list_of_strings is None:
        return list_of_strings

    list_of_floats = [float(string) for string in list_of_strings]
    return list_of_floats


def cast_list_to_bool(
        list_of_strings: Union[List[str], None]) -> Union[List[bool], None]:
    """Converts a list of strings into a list of booleans.

    Args:
        list_of_strings: List containing string variables.

    Returns:
        List of booleans or None, if the list of strings is None.
    """
    if list_of_strings is None:
        return list_of_strings

    for string in list_of_strings:
        if string not in ["True", "False", True, False]:
            raise ValueError(f"Could not cast string `{string}` to boolean.")

    list_of_bools = [(str(string) == "True") or (str(string) != "False")
                     for string in list_of_strings]
    return list_of_bools
