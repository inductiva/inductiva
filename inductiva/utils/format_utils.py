"""Util functions for formatting data for printing to console."""
from distutils.util import strtobool
import datetime
import os

import pandas as pd
import numpy as np


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


def get_tabular_str(
    rows,
    columns,
    default_col_space=15,
    override_col_space=None,
    formatters=None,
) -> str:
    """Converts a list of lists to a string table.

    Temporary solution to display tables. `pandas` dependency to print tables
    is overkill, but since we have `pandas` anyway for now, we can use it.
    """
    df = pd.DataFrame(rows, columns=columns)
    formatters = formatters or {}

    col_space = {col: default_col_space for col in columns}
    col_space.update(override_col_space or {})

    # replace None with np.nan so that pandas can format them as "n/a"
    # by passing na_rep="n/a" to to_string()
    df.fillna(np.nan, inplace=True)

    return df.to_string(
        index=False,
        na_rep="n/a",
        formatters=formatters,
        col_space=col_space,
    )


def get_dataframe_str(dataframe,
                      default_col_space=15,
                      override_col_space=None,
                      formatters=None):
    """Converts a dataframe to a string table.

    Temporary solution to display tables. `pandas` dependency to print tables
    is overkill, but since we have `pandas` anyway for now, we can use it.
    """
    formatters = formatters or {}

    col_space = {col: default_col_space for col in dataframe.columns}
    col_space.update(override_col_space or {})

    # replace None with np.nan so that pandas can format them as "n/a"
    # by passing na_rep="n/a" to to_string()
    dataframe.fillna(np.nan, inplace=True)

    return dataframe.to_string(index=False,
                               na_rep="n/a",
                               formatters=formatters,
                               col_space=col_space)
