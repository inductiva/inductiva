"""Register CLI command for logout."""
import argparse
import os

from inductiva import constants


def logout(_):
    """
    Logs the user out by removing the stored API key.

    The function will remove the stored API key, logging the user out.
    """
    if not os.path.exists(constants.API_KEY_FILE_PATH):
        print("Error: You are not logged in.")
        return

    os.remove(constants.API_KEY_FILE_PATH)
    print("Logout successful.")


def register(parser):
    """Register the logout command."""
    subparser = parser.add_parser(
        "logout",
        help="Logout by removing the Inductiva API key.",
        formatter_class=argparse.RawTextHelpFormatter)

    subparser.description = (
        "The `inductiva logout` command allows you to log out by removing "
        "the stored API key.\n")

    subparser.set_defaults(func=logout)
