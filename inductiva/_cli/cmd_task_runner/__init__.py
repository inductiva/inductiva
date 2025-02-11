"""Register CLI commands for storage."""
import argparse
import os

from inductiva import constants
from inductiva._cli import loader, utils


def register(root_parser):

    parser = root_parser.add_parser(
        "task-runner",
        help="Task-Runner management utilities.",
        formatter_class=argparse.RawTextHelpFormatter)

    parser.description = (
        "Task-Runner management utilities.\n\n"
        "The `inductiva task-runner` command allows you to "
        "manage your local task-runners on the platform.\n"
        "It provides utilities for launching or terminating task-runners.\n")

    utils.show_help_msg(parser)

    subparsers = parser.add_subparsers(title="available subcomands")
    loader.load_commands(subparsers,
                         os.path.dirname(__file__),
                         package=__name__,
                         ignores_prefix=constants.LOADER_IGNORE_PREFIX,
                         hides_prefix=constants.LOADER_HIDE_PREFIX)
