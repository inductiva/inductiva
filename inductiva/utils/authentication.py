"""This module contains utility functions for handling authentication. """

import re

from inductiva import constants


def is_valid_token(token):
    """Validate token characters."""
    pattern = re.compile(r"^[A-Za-z0-9._-]+$")
    if not pattern.match(token):
        return False
    return True


def get_stored_api_key(path=None):
    """Returns the stored API key from the file system."""
    key_path = path or constants.API_KEY_FILE_PATH
    if not key_path.exists():
        return ""
    with open(key_path, "r", encoding="utf-8") as f:
        return f.read()


def set_stored_api_key(api_key, path=None) -> bool:
    """Saves the API key to the file system.
    Returns True if the key already existed, False otherwise.
    """
    key_path = path or constants.API_KEY_FILE_PATH

    exists = key_path.exists()

    with open(key_path, "w", encoding="utf-8") as f:
        f.write(api_key)
    return exists


def remove_stored_api_key(path=None):
    """Removes the stored API key from the file system."""
    key_path = path or constants.API_KEY_FILE_PATH
    if key_path.exists():
        key_path.unlink()
