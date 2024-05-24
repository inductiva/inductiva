"""Register CLI commands for projetcs."""
import argparse
import os

from inductiva import constants
from inductiva._cli import loader, utils


def register(root_parser):

    parser = root_parser.add_parser(
        "projects",
        help="Projects management utilities.",
        formatter_class=argparse.RawTextHelpFormatter)

    parser.description = ("Projects management utilities.\n\n"
                          "The `inductiva projects` command allows you to "
                          "consult existing projetcs.\n")

    utils.show_help_msg(parser)

    subparsers = parser.add_subparsers(title="available subcomands")
    loader.load_commands(subparsers,
                         os.path.dirname(__file__),
                         package=__name__,
                         ignores_prefix=constants.LOADER_IGNORE_PREFIX)
