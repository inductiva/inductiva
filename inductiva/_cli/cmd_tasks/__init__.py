"""Register CLI commands for storage."""
import argparse
import os
import textwrap

from inductiva import constants
from inductiva._cli import loader, utils


def register(root_parser):

    parser = root_parser.add_parser(
        "tasks",
        help="Task management utilities.",
        formatter_class=argparse.RawTextHelpFormatter)

    parser.description = textwrap.dedent("""\
        Task management utilities.

        The `inductiva tasks` command provides tools to manage your tasks on 
        the Inductiva platform. It includes utilities for monitoring task 
        progress, viewing details, and terminating tasks when needed.
    """)

    utils.show_help_msg(parser)

    subparsers = parser.add_subparsers(title="available subcomands")
    loader.load_commands(subparsers,
                         os.path.dirname(__file__),
                         package=__name__,
                         ignores_prefix=constants.LOADER_IGNORE_PREFIX,
                         hides_prefix=constants.LOADER_HIDE_PREFIX)
