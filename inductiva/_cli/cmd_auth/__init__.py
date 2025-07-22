"""Register CLI commands for logs."""
import argparse
import os
import textwrap

from .. import loader, utils
from ... import constants


def register(root_parser):
    parser = root_parser.add_parser(
        "auth",
        help="Authentication commands",
        formatter_class=argparse.RawTextHelpFormatter)

    parser.description = textwrap.dedent("""\
        Authentication management utilities.

        The `inductiva auth` command allows you to manage your authentication
        on Inductiva API.

        This is a fundamental aspect of your experience: you will only be able 
        to start machines and launch simulations after you are authenticated. 
        Once authenticated, your credentials will be stored locally for future 
        sessions.
                                         
        However, you will need to authenticate separately on each local machine
        from which you want to use Inductiva.
    """)

    utils.show_help_msg(parser)

    subparsers = parser.add_subparsers(title="available subcomands")
    loader.load_commands(subparsers,
                         os.path.dirname(__file__),
                         package=__name__,
                         ignores_prefix=constants.LOADER_IGNORE_PREFIX)
