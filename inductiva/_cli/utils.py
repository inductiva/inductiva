"""CLI utils."""
import os

import inductiva
from inductiva import constants


def check_running_for_first_time():
    """Checks if it is the first time the cli is being called.

    It does this by checking if there is a folder:

      `$HOME/.inductiva/v{inductiva_version}`

    If there is not returns true and creates that folder.

    """
    version = inductiva.__version__
    dir_name = constants.LOCAL_LOGGING_DIR / f"v{version}"
    if not os.path.exists(dir_name):
        os.mkdir(dir_name)
        return True
    return False


def show_help_msg(parser):
    """Show help message for command if no subcommand is provided."""
    parser.set_defaults(func=lambda _: parser.print_help())
