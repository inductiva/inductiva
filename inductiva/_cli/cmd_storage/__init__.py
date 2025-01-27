"""Register CLI commands for storage."""
import argparse
import os

from inductiva import constants
from inductiva._cli import loader, utils


def register(root_parser):

    parser = root_parser.add_parser(
        "storage",
        help="Remote storage management utilities.",
        formatter_class=argparse.RawTextHelpFormatter)
    parser.description = (
        "Remote storage management utilities.\n\n"
        "The `inductiva storage` command oversees your data on the platform.\n"
        "It enables listing all items, removing specific ones,\n"
        "and calculating the total size of your storage.\n")
    utils.show_help_msg(parser)

    subparsers = parser.add_subparsers(title="available subcomands")
    loader.load_commands(subparsers,
                         os.path.dirname(__file__),
                         package=__name__,
                         ignores_prefix=constants.LOADER_IGNORE_PREFIX,
                         hides_prefix=constants.LOADER_HIDE_PREFIX)
