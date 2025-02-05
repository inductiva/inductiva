"""Register CLI commands for projects."""
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
                          "The `inductiva projects list` command allows you to "
                          "consult existing projects.\n"
                          "The `inductiva projects download` command allows you"
                          " to download the files of the tasks of a project.\n")

    utils.show_help_msg(parser)

    subparsers = parser.add_subparsers(title="available subcomands")
    loader.load_commands(subparsers,
                         os.path.dirname(__file__),
                         package=__name__,
                         ignores_prefix=constants.LOADER_IGNORE_PREFIX,
                         hides_prefix=constants.LOADER_HIDE_PREFIX)
