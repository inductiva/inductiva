"""Register CLI commands for autocomplete."""
import argparse
import os

from inductiva import constants
from inductiva._cli import loader, utils


def register(root_parser):
    parser = root_parser.add_parser(
        "autocomplete",
        help="Enables auto-completion for inductiva commands.",
        formatter_class=argparse.RawTextHelpFormatter)

    parser.description = (
        "Controls the autocomplete behavior for the inductiva package.\n"
        "The `inductiva autocomplete` command allows you to setup "
        "autocompletion for the cli.\n"
        "At the moment we only provide support for zsh.")

    utils.show_help_msg(parser)

    subparsers = parser.add_subparsers(title="available subcomands")
    loader.load_commands(subparsers,
                         os.path.dirname(__file__),
                         package=__name__,
                         ignores_prefix=constants.LOADER_IGNORE_PREFIX)
