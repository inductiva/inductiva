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


#clean this function
def apply_formatters(rows, columns, formatters):
    data = {}

    for index, column_name in enumerate(columns):
        data[column_name] = [row[index] for row in rows]

    for column_name, formatter in formatters.items():
        if column_name in data:
            data[column_name] = [formatter(x) for x in data[column_name]]

    return data


def get_tabular_str(
    rows,
    columns,
    formatters=None,
) -> str:
    """Converts a list of lists to a string table.

    """

    formatters = formatters or {}

    data = apply_formatters(rows, columns, formatters)
    data_tabulated = tabulate(data, headers=columns, missingval="n/a")

    return data_tabulated


def get_tasks_str(columns, rows, formatters=None):
    """Converts a list of tasks to a list of ids.

    """

    formatters = formatters or {}

    data = apply_formatters(rows, columns, formatters)

    data_tabulated = tabulate(data, headers=columns, missingval="n/a")

    # replace None with np.nan so that pandas can format them as "n/a"
    # by passing na_rep="n/a" to to_string()

    return data_tabulated
