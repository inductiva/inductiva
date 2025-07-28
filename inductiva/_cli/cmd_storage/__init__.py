"""Register CLI commands for storage."""
import argparse
import os
import textwrap

from inductiva import constants
from inductiva._cli import loader, utils


def register(root_parser):

    parser = root_parser.add_parser(
        "storage",
        help="Remote storage management utilities.",
        formatter_class=argparse.RawTextHelpFormatter)

    parser.description = textwrap.dedent("""\
        Remote storage management utilities. 
                                        
        Use the `inductiva storage` command to manage your data in Inductiva's
        remote storage. It supports a range of operations including listing
        files, downloading and deleting items, and calculating total storage
        usage and cost.
    """)

    utils.show_help_msg(parser)

    subparsers = parser.add_subparsers(title="available subcomands")
    loader.load_commands(subparsers,
                         os.path.dirname(__file__),
                         package=__name__,
                         ignores_prefix=constants.LOADER_IGNORE_PREFIX,
                         hides_prefix=constants.LOADER_HIDE_PREFIX)
