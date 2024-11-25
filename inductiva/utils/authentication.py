"""This module contains utility functions for handling authentication. """

from inductiva import constants


def get_stored_api_key(path=None):
    """Returns the stored API key from the file system."""
    key_path = path or constants.API_KEY_FILE_PATH
    if not key_path.exists():
        return ""
    with open(key_path, "r", encoding="utf-8") as f:
        return f.read()


def set_stored_api_key(api_key, path=None):
    """Saves the API key to the file system."""
    key_path = path or constants.API_KEY_FILE_PATH
    with open(key_path, "w", encoding="utf-8") as f:
        f.write(api_key)


def remove_stored_api_key(path=None):
    """Removes the stored API key from the file system."""
    key_path = path or constants.API_KEY_FILE_PATH
    if key_path.exists():
        key_path.unlink()
