"""Util functions for formatting data for printing to console."""
import datetime
import pandas as pd
import numpy as np


def bytes_formatter(n_bytes: int) -> str:
    """Convert bytes to human readable string."""
    res = float(n_bytes)

    for unit in ["B", "KiB", "MiB", "GiB", "TiB"]:
        if res < 1024:
            if unit == "B":
                return f"{res:.0f} {unit}"
            else:
                return f"{res:.2f} {unit}"
        res /= 1024

    return f"{bytes_formatter:.2f} PiB"


def datetime_formatter(dt: str) -> str:
    # get time in local timezone
    local_dt = datetime.datetime.fromisoformat(dt).astimezone()
    return local_dt.strftime("%d %b, %H:%M:%S")


def seconds_formatter(secs: float) -> str:
    return (f"{int(secs // 3600)}h "
            f"{int((secs % 3600) // 60)}m "
            f"{int(secs % 60)}s")


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
