"""CLI utils."""
import argparse
import os

import inductiva
from inductiva import constants

BEGIN_TAG = "# >>> INDUCTIVA BEGIN:"
END_TAG = "# >>> INDUCTIVA END:"


def positive_float(value):
    """Check if the given value is a positive float."""
    try:
        value = float(value)
    except TypeError as _:
        raise argparse.ArgumentTypeError(f"{value} is an invalid type.")
    if value <= 0:
        raise argparse.ArgumentTypeError(f"{value} must be a positive number.")
    return value


def add_watch_argument(subparser):
    """Add watch flag to the subparser that prompts the commands every Nsecs."""
    subparser.set_defaults(watchable=True)
    subparser.add_argument("-w",
                           "--watch",
                           nargs="?",
                           const=2.0,
                           type=positive_float,
                           help="Prompt the command every N seconds.")


def check_running_for_first_time():
    """Checks if it is the first time the cli is being called.

    It does this by checking if there is a folder:

      `$HOME/.inductiva/v{inductiva_version}`

    If there is not returns true and creates that folder.

    """
    version = inductiva.__version__
    dir_name = constants.HOME_DIR / f"v{version}"
    if not os.path.exists(dir_name):
        os.mkdir(dir_name)
        return True
    return False


def show_help_msg(parser):
    """Show help message for command if no subcommand is provided."""
    parser.set_defaults(func=lambda _: parser.print_help())
