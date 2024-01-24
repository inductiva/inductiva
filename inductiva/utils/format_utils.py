"""Util functions for formatting data for printing to console."""
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
    """Apply formatters to a list of lists of data (rows).

    Args:
        table_data (dict): Dictionary of column names and lists of data.
        formatters (dict): Dictionary of column names and functions to
            apply to that column's data.
    """
    for column_name, formatter in formatters.items():
        if column_name in table_data:
            table_data[column_name] = [
                formatter(x) for x in table_data[column_name]
            ]

    return table_data


def get_tabular_str(rows: list, columns: list, formatters: dict = None) -> str:
    """Converts the list rows to a string table.

    Args:
        rows (list): A list of data to be printed.
        columns (list): A list of column names.
        formatters (dict): A dictionary of column names and functions to apply
            to the data in that column. The function should take a single
            argument and return a string. The function will be applied to the
            data in the column before printing. Defaults to None.

    """

    formatters = formatters or {}

    data = {}

    for index, column_name in enumerate(columns):
        data[column_name] = [row[index] for row in rows]

    data = apply_formatters(data, formatters)
    data_tabulated_str = tabulate(data, headers=columns, missingval="n/a")

    return data_tabulated_str
