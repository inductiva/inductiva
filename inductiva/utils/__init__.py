# pylint: disable=missing-module-docstring
from . import files
from . import data
from .files import download_from_url
from .authentication import (
    get_stored_api_key,
    is_valid_token,
    remove_stored_api_key,
    set_stored_api_key,
)


def strtobool(val: str) -> bool:
    """Convert a string representation of truth to a boolean."""
    if isinstance(val, bool):
        return val

    if isinstance(val, str):
        val = val.lower()
        if val in ("y", "yes", "t", "true", "on", "1"):
            return True
        elif val in ("n", "no", "f", "false", "off", "0"):
            return False

    raise ValueError(f"Invalid truth value: {val}")
