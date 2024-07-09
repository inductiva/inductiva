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

    print(f"■ Name: {user_info['name']}", file=fout)
    print(f"■ Email: {user_info['email']}", file=fout)
    print(f"■ Username: {user_info['username']}", file=fout)

    tier = user_info["tier"]["name"]
    available_credits = user_info["tier"]["available_credits"]
    print(f"     ■ Tier: {tier}", file=fout)
    print(f"     ■ Available credits: {available_credits}", file=fout)

    for program in user_info["programs"]:
        table["name"].append(program["name"])
        table["enrollment date"].append(program["enrollment_date"])
        table["expiry date"].append(program["expiry_date"])
        table["available credits"].append(program["available_credits"])
        table["initial credits"].append(program["initial_credits"])

    emph_formatter = format_utils.get_ansi_formatter()

    header_formatters = [
        lambda x: emph_formatter(x.upper(), format_utils.Emphasis.BOLD)
    ]

    formatters = {
        "enrollment date": [format_utils.datetime_formatter,],
        "expiry date": [format_utils.datetime_formatter,]
    }

    table = format_utils.get_tabular_str(table,
                                         formatters=formatters,
                                         header_formatters=header_formatters)
    print("════════════════════════════════════", file=fout)
    print("Programs:", file=fout)
    print(table, file=fout)

    total_available_credits = user_info["total_available_credits"]
    print("════════════════════════════════════", file=fout)
    print(f"Total available credits: {total_available_credits}", file=fout)
    print("════════════════════════════════════", file=fout)
    return 0


def register(parser):
    """Register the user credits command."""

    subparser = parser.add_parser("info",
                                  help="List the user's information.",
                                  formatter_class=argparse.RawTextHelpFormatter)

    subparser.description = ("The `inductiva user info` command provides "
                             "an overview of your tier, programs and credits.\n"
                             "It lists all your programs as well as the"
                             "credits left for you to use.\n")

    _cli.utils.add_watch_argument(subparser)

    subparser.set_defaults(func=get_credits)
