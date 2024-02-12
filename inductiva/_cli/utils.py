"""CLI utils."""
import os

import inductiva
from inductiva import constants


def check_running_for_first_time():
    version = inductiva.__version__.replace(".", "-")
    dir_name = constants.LOCAL_LOGGING_DIR / f"v{version}"
    if not os.path.exists(dir_name):
        os.mkdir(dir_name)
        return True
    return False


def show_help_msg(parser):
    """Show help message for command if no subcommand is provided."""
    parser.set_defaults(func=lambda _: parser.print_help())
