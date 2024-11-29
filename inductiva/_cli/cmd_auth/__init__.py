"""Register CLI commands for logs."""
import argparse
import os

from .. import loader, utils
from ... import constants


def register(root_parser):
    parser = root_parser.add_parser(
        "auth",
        help="Authentication commands",
        formatter_class=argparse.RawTextHelpFormatter)

    parser.description = ("Authentication management utilities.\n\n"
                          "The `inductiva auth` command allows you to "
                          "manage user authentication.\n")

    utils.show_help_msg(parser)

    subparsers = parser.add_subparsers(title="available subcomands")
    loader.load_commands(subparsers,
                         os.path.dirname(__file__),
                         package=__name__,
                         ignores_prefix=constants.LOADER_IGNORE_PREFIX)
