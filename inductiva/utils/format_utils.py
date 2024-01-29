"""Util functions for formatting data for printing to console."""
from typing import Any, Callable, Dict, Iterable, Mapping, Optional, Union
from distutils.util import strtobool
import datetime
import os

from tabulate import tabulate


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
    for column_name, formatter in formatters.items():
        if column_name in table_data:
            table_data[column_name] = [
                formatter(x) for x in table_data[column_name]
            ]

    return table_data


def get_tabular_str(tabular_data: Union[Mapping[str, Iterable[Any]],
                                        Iterable[Iterable[Any]]],
                    headers: Optional[Iterable[Any]] = None,
                    formatters: Optional[Dict[str, Callable]] = None) -> str:
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
        headers = tabular_data.keys()

    tabular_data_formatted = apply_formatters(tabular_data, formatters)

    return tabulate(tabular_data_formatted, headers=headers, missingval="n/a")
