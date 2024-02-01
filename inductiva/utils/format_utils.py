"""Util functions for formatting data for printing to console."""
from typing import (Any, Callable, Dict, Iterable, Mapping, Optional, Tuple,
                    Union, List)
from distutils.util import strtobool
import datetime
import os
import sys
import copy

import tabulate

tabulate.PRESERVE_WHITESPACE = True

EMPHASIS = {
    "red": ("\033[31m", "\033[0m"),
    "green": ("\033[92m", "\033[0m"),
    "bold": ("\u001b[1m", "\u001b[0m")
}


def getenv_bool(varname, default):
    """Get boolean value from environment variable."""
    return bool(strtobool(os.getenv(varname, str(default))))


def bytes_formatter(n_bytes: int) -> str:
    """Convert bytes to human readable string."""
    res = float(n_bytes)

    for unit in ["B", "KB", "MB", "GB", "TB"]:
        if res < 1000:
            if unit == "B":
                return f"{res:.0f} {unit}"
            else:
                return f"{res:.2f} {unit}"
        res /= 1000

    return f"{bytes_formatter:.2f} PB"


def _supports_ansi():
    if sys.platform.startswith("win"):
        return "TERM" in os.environ and os.environ["TERM"] == "xterm"
    return hasattr(sys.stdout, "isatty") and sys.stdout.isatty()


def emphasis_formater(string_to_emphasize, emphasis):
    if not _supports_ansi():
        return string_to_emphasize
    if not emphasis in ["red", "green", "bold"]:
        raise ValueError("Desired emphasis is not supported. "
                         "Select either `red`, `green` or `bold`")
    emph = EMPHASIS[emphasis]
    return f"{emph[0]}{string_to_emphasize}{emph[1]}"


def spacing_formater(x, num_spaces=6):
    return f"{x}" + num_spaces * " "


def datetime_formatter(dt: str) -> str:
    # get time in local timezone
    if dt is None:
        return None
    local_dt = datetime.datetime.fromisoformat(dt).astimezone()
    return local_dt.strftime("%d %b, %H:%M:%S")


def seconds_formatter(secs: float) -> str:
    """Convert seconds to time human readable string."""
    return str(datetime.timedelta(seconds=round(secs)))


def apply_formatters(table_data: dict, formatters: dict):
    """Applies a dict of formatters to dict of data.

    Args:
        table_data : Dictionary of column names and lists of data.
        formatters : Dictionary of column names and functions to
            apply to that column's data.
    """
    output_table_data = copy.deepcopy(table_data)
    for column_name, formatters in formatters.items():
        if column_name in output_table_data:
            for formatter in formatters:
                output_table_data[column_name] = [
                    formatter(x) for x in output_table_data[column_name]
                ]

    return output_table_data


def get_tabular_data(
    tabular_data: Union[Mapping[str, Iterable[Any]], Iterable[Iterable[Any]]],
    headers: Optional[Iterable[Any]] = None,
    formatters: Optional[Dict[str,
                              List[Callable]]] = None) -> Tuple[dict, list]:
    """Converts a table of data (Mapping or any Iterable) to
    dict and a list of headers.

    Args:
        Gets the same arguments as get_tabular_str.
    Returns:
        A dict of column names and iterable of data, and a list of headers.
    """
    formatters = formatters or {}
    headers = headers or []

    if not isinstance(tabular_data, Mapping):

        #if we have no headers data will be empty.
        #So, we want our original tabular_data
        if headers:
            tabular_data = {
                header: [row[index] for row in tabular_data]
                for index, header in enumerate(headers)
            }
    else:
        headers = list(tabular_data.keys())

    tabular_data_formatted = apply_formatters(tabular_data, formatters)
    return tabular_data_formatted, headers


def _table_indenter(table_string, num_spaces):
    "Adds spaces to the beggining of each table line."
    indentenation = " " * num_spaces
    return ("\n" + indentenation).join(table_string.split("\n"))


def get_tabular_str(tabular_data: Union[Mapping[str, Iterable[Any]],
                                        Iterable[Iterable[Any]]],
                    headers: Optional[Iterable[Any]] = None,
                    formatters: Optional[Dict[str, Callable]] = None,
                    header_formatters: Optional[List[Callable]] = None,
                    indentation_level: Optional[int] = 4) -> str:
    """Converts a table of data (Mapping or any Iterable) to a string table.

    Args:
        tabular_data: can be a list-of-lists (or another iterable of 
            iterables), a list of named tuples, a dictionary of
            iterables, an iterable of dictionaries, an iterable of
            dataclasses (Python 3.7+)
        headers: A list of column names.
            Only needed if tabular_data is not a Mapping.
        formatters: A dictionary of column names and functions to apply
            to the data in that column. The function should take a single
            argument and return a string. The function will be applied to the
            data in the column before printing. Defaults to None.
    Returns:
        A string table with the contents of tabular_data and headers.
    """

    data, headers = get_tabular_data(tabular_data, headers, formatters)

    for formatter in header_formatters:
        headers = [formatter(header) for header in headers]
    if indentation_level is not None:
        headers = [" " * indentation_level + header for header in headers]

    table = tabulate.tabulate(data,
                              headers=headers,
                              missingval="n/a",
                              tablefmt="plain")
    if indentation_level is not None:
        table = _table_indenter(table, indentation_level)

    return f"\n{table}\n"
