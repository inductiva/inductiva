"""Register CLI commands for simulators."""
import argparse
import os

from inductiva._cli import loader, utils
from inductiva import constants


def register(root_parser):
    """Register the "simulators" sub-command."""

    parser = root_parser.add_parser(
        "simulators",
        help="Information about available simulators.",
        formatter_class=argparse.RawTextHelpFormatter)

    parser.description = (
        "Information about available simulators.\n\n"
        "The `inductiva simulators` command provides utility sub-commands\n"
        "for managing the available simulators within the Inductiva API.\n")

    utils.show_help_msg(parser)
    subparsers = parser.add_subparsers(title="available subcomands")
    loader.load_commands(subparsers,
                         os.path.dirname(__file__),
                         package=__name__,
                         ignores_prefix=constants.LOADER_IGNORE_PREFIX,
                         hides_prefix=constants.LOADER_HIDE_PREFIX)
