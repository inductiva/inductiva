"""This module contains utility functions for handling authentication. """

import base64
import json
import re

from inductiva import constants


def is_valid_jwe(token):
    try:
        # Split the token into its parts
        parts = token.split(".")
        if len(parts) != 5:
            return False

        # Check if the token parts contain only valid base64url characters
        base64url_pattern = re.compile(r'^[A-Za-z0-9_-]+$')
        if not all(base64url_pattern.match(part) for part in parts):
            return False

        # Decode the header
        header = base64.urlsafe_b64decode(parts[0] + '==')
        header_json = json.loads(header)

        # Check if the header contains the required fields
        if "alg" not in header_json or "enc" not in header_json:
            return False

        return True
    except (ValueError, json.JSONDecodeError, base64.binascii.Error):
        return False


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
