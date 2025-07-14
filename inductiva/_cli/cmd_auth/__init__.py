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

    parser.description = \
"""
Authentication management utilities.

The `inductiva auth` command allows you to manage user authentication on
Inductiva API.

This is a fundamental aspect of your experience: you will only be able to
start machines and launch simulations after you are authenticated. Once
authenticated, your credentials will be stored locally for future sessions.
Howecer, you will need to perform the authentication step from every local
machine you want to use Inductiva from.
"""

    utils.show_help_msg(parser)

    subparsers = parser.add_subparsers(title="available subcomands")
    loader.load_commands(subparsers,
                         os.path.dirname(__file__),
                         package=__name__,
                         ignores_prefix=constants.LOADER_IGNORE_PREFIX)
