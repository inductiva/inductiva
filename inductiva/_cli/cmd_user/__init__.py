"""Register CLI commands for user."""
import argparse
import os
import textwrap

from inductiva._cli import loader, utils
from inductiva import constants


def register(root_parser):
    parser = root_parser.add_parser(
        "user",
        help="User management utilities.",
        formatter_class=argparse.RawTextHelpFormatter)

    parser.description = textwrap.dedent("""\
        User management utilities.

        The `inductiva user` command allows you to consult information
        about your Inductiva account, including available credits and active
        quotas.                              
    """)

    utils.show_help_msg(parser)

    subparsers = parser.add_subparsers(title="available subcomands")
    loader.load_commands(subparsers,
                         os.path.dirname(__file__),
                         package=__name__,
                         ignores_prefix=constants.LOADER_IGNORE_PREFIX,
                         hides_prefix=constants.LOADER_HIDE_PREFIX)
