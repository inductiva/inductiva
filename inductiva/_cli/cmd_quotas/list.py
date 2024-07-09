"""List the user quotas information via CLI."""
from typing import TextIO
import argparse
import logging
import sys

from inductiva import _cli
from inductiva._cli.cmd_user import info


def get_quotas(_, fout: TextIO = sys.stdout):
    """ Lists the user's quotas.

    Lists all the user's quotas and the quotas left for the user to use.
    """
    logging.warning("\nWARNING: `inductiva quotas list` is deprecated"
                    " and will be removed in a future release.\n"
                    "Please use the `inductiva user info` command instead.\n")
    info.get_info(_, fout)
    return 0


def register(parser):
    """Register the quotas list command."""
    subparser = parser.add_parser("list",
                                  help="List the user's quotas.",
                                  formatter_class=argparse.RawTextHelpFormatter)

    subparser.description = ("The `inductiva quotas list` command provides "
                             "an overview of your quotas.\n"
                             "It lists all your available quotas as well as the"
                             "quotas left for you to use\n")

    _cli.utils.add_watch_argument(subparser)

    subparser.set_defaults(func=get_quotas)
