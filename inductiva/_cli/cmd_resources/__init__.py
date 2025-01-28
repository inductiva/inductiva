"""Register CLI commands for logs."""
import os
import argparse

from inductiva import constants
from inductiva._cli import loader, utils


def register(root_parser):

    parser = root_parser.add_parser(
        "resources",
        help="Computational resource management utilities.",
        formatter_class=argparse.RawTextHelpFormatter)
    parser.description = (
        "Computational resource management utilities.\n\n"
        "The `inductiva resources` command provides utilities "
        "for managing computational resources.\nIt allows you "
        "to show estimated costs of resources, "
        "show available machine types, list current resources\n"
        " being used, and terminate resources.\n")

    utils.show_help_msg(parser)

    subparsers = parser.add_subparsers(title="available subcomands")
    loader.load_commands(subparsers,
                         os.path.dirname(__file__),
                         package=__name__,
                         ignores_prefix=constants.LOADER_IGNORE_PREFIX,
                         hides_prefix=constants.LOADER_HIDE_PREFIX)
