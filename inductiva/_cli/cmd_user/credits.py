"""List the user credits information via CLI."""
from collections import defaultdict
from typing import TextIO
import argparse
import sys

from inductiva import users, _cli
from inductiva.utils import format_utils


def get_credits(_, fout: TextIO = sys.stdout):
    """ Lists the user's credits.

    Lists all the user's credits and the credits left for the user to use.
    """
    table = defaultdict(list)

    user_info = users.get_info()
    tier = user_info["tier"]
    print(f"Tier: {tier}", file=fout)

    for name, credit in user_info["credits"].items():
        table[""].append(name)
        table["current usage"].append(credit["used"])
        table["max allowed"].append(credit["total"])
        table["remaining"].append(credit["remaining"])

    emph_formatter = format_utils.get_ansi_formatter()

    header_formatters = [
        lambda x: emph_formatter(x.upper(), format_utils.Emphasis.BOLD)
    ]

    table = format_utils.get_tabular_str(table,
                                         header_formatters=header_formatters)

    print(table, file=fout, end="")

    return 0


def register(parser):
    """Register the user credits command."""

    subparser = parser.add_parser("credits",
                                  help="List the user's credits.",
                                  formatter_class=argparse.RawTextHelpFormatter)

    subparser.description = ("The `inductiva user credits` command provides "
                             "an overview of your tier and credits.\n"
                             "It lists all your total credits as well as the"
                             "credits left for you to use.\n")

    _cli.utils.add_watch_argument(subparser)

    subparser.set_defaults(func=get_credits)
