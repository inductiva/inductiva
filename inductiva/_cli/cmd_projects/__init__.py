"""Register CLI commands for projects."""
import argparse
import os
import textwrap

from inductiva import constants
from inductiva._cli import loader, utils


def register(root_parser):

    parser = root_parser.add_parser(
        "projects",
        help="Projects management utilities.",
        formatter_class=argparse.RawTextHelpFormatter)

    parser.description = textwrap.dedent("""\
        Projects management utilities.

        The `inductiva projects` command offers tools to manage your projects 
        and access related information, such as listing all projects and 
        downloading their associated task output files.
    """)

    utils.show_help_msg(parser)

    subparsers = parser.add_subparsers(title="available subcomands")
    loader.load_commands(subparsers,
                         os.path.dirname(__file__),
                         package=__name__,
                         ignores_prefix=constants.LOADER_IGNORE_PREFIX,
                         hides_prefix=constants.LOADER_HIDE_PREFIX)
