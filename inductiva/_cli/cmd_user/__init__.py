"""Register CLI commands for user."""
import argparse
import os

from inductiva._cli import loader, utils
from inductiva import constants


def register(root_parser):
    parser = root_parser.add_parser(
        "user",
        help="User management utilities.",
        formatter_class=argparse.RawTextHelpFormatter)

    parser.description = ("User management utilities.\n\n"
                          "The `inductiva user` command allows you to "
                          "consult your user's internal information.\n")

    utils.show_help_msg(parser)

    subparsers = parser.add_subparsers(title="available subcomands")
    loader.load_commands(subparsers,
                         os.path.dirname(__file__),
                         package=__name__,
                         ignores_prefix=constants.LOADER_IGNORE_PREFIX,
                         hides_prefix=constants.LOADER_HIDE_PREFIX)
