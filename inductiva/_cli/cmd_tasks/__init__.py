"""Register CLI commands for storage."""
import argparse
import os

from inductiva import constants
from inductiva._cli import loader, utils


def register(root_parser):

    parser = root_parser.add_parser(
        "tasks",
        help="Task management utilities.",
        formatter_class=argparse.RawTextHelpFormatter)

    parser.description = (
        "Task management utilities.\n\n"
        "The `inductiva tasks` command allows you to "
        "manage your tasks on the platform.\n"
        "It provides utilities for monitoring and terminating tasks.\n")

    utils.show_help_msg(parser)

    subparsers = parser.add_subparsers(title="available subcomands")
    loader.load_commands(subparsers,
                         os.path.dirname(__file__),
                         package=__name__,
                         ignores_prefix=constants.LOADER_IGNORE_PREFIX,
                         hides_prefix=constants.LOADER_HIDE_PREFIX)
