"""Register CLI commands for storage."""
from inductiva import constants
import os

from inductiva._cli import loader, utils


def register(root_parser):

    parser = root_parser.add_parser("tasks", help="Task management utilities.")
    utils.show_help_msg(parser)

    subparsers = parser.add_subparsers()
    loader.load_commands(subparsers,
                         os.path.dirname(__file__),
                         package=__name__,
                         ignores_prefix=constants.LOADER_IGNORE_PREFIX)
