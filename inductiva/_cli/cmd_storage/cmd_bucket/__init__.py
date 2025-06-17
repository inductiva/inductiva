"""Register CLI commands for storage buckets."""
import argparse
import os

from inductiva import constants
from inductiva._cli import loader, utils


def register(root_parser):
    parser = root_parser.add_parser(
        "bucket",
        help="Storage bucket management utilities.",
        formatter_class=argparse.RawTextHelpFormatter,
    )

    parser.description = (
        "Storage bucket management utilities.\n\n"
        "The `inductiva storage bucket` command provides tools for "
        "managing with cloud storage buckets.\n"
        "Use it to manage cloud storage buckets linked to your account.")

    utils.show_help_msg(parser)

    subparsers = parser.add_subparsers(title="available subcomands")
    loader.load_commands(subparsers,
                         os.path.dirname(__file__),
                         package=__name__,
                         ignores_prefix=constants.LOADER_IGNORE_PREFIX,
                         hides_prefix=constants.LOADER_HIDE_PREFIX)
