"""Register CLI commands for quotas."""
import argparse
import os

from inductiva import constants
from inductiva._cli import loader, utils


def register(root_parser):
    parser = root_parser.add_parser(
        "quotas",
        help="Quotas management utilities.",
        formatter_class=argparse.RawTextHelpFormatter)

    parser.description = ("Quotas management utilities.\n\n"
                          "The `inductiva quotas` command allows you to "
                          "consult your user's internal quotas.\n")

    utils.show_help_msg(parser)

    subparsers = parser.add_subparsers(title="available subcomands")
    loader.load_commands(subparsers,
                         os.path.dirname(__file__),
                         package=__name__,
                         ignores_prefix=constants.LOADER_IGNORE_PREFIX,
                         hides_prefix=constants.LOADER_HIDE_PREFIX)
