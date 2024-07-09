"""List the user quotas information via CLI."""
from collections import defaultdict
from typing import TextIO
import argparse
import sys

from inductiva import users, _cli
from inductiva.utils import format_utils


def get_quotas(_, fout: TextIO = sys.stdout):
    """ Lists the user's quotas.

    Lists all the user's quotas and the quotas left for the user to use.
    """
    table = defaultdict(list)

    table_instance = defaultdict(list)

    for _, quota in users.get_quotas().items():
        current = quota['in_use']
        max_allowed = quota['max_allowed']

        current_str = (f"{current} {quota['unit']}"
                       if current is not None else "N/A")
        max_allowed_str = (f"{max_allowed} {quota['unit']}"
                           if max_allowed is not None else "N/A")

        if quota["scope"] == "instance":
            table_instance[""].append(quota["label"])
            table_instance["max allowed"].append(max_allowed_str)
        else:
            table[""].append(quota["label"])
            table["current usage"].append(current_str)
            table["max allowed"].append(max_allowed_str)

    emph_formatter = format_utils.get_ansi_formatter()

    header_formatters = [
        lambda x: emph_formatter(x.upper(), format_utils.Emphasis.BOLD)
    ]

    table = format_utils.get_tabular_str(table,
                                         header_formatters=header_formatters)

    table_instance = format_utils.get_tabular_str(
        table_instance, header_formatters=header_formatters)
    print("■ Global User quotas", end="")
    print(table, file=fout)
    print("═══════════════════════════════════════════════════════════════════")
    print("■ Instance User quotas", end="")
    print(table_instance, file=fout, end="")

    return 0


def register(parser):
    """Register the user quotas command."""
    subparser = parser.add_parser("quotas",
                                  help="List the user's quotas.",
                                  formatter_class=argparse.RawTextHelpFormatter)

    subparser.description = ("The `inductiva user quotas` command provides "
                             "an overview of your quotas.\n"
                             "It lists all your available quotas as well as the"
                             "quotas left for you to use\n")

    _cli.utils.add_watch_argument(subparser)

    subparser.set_defaults(func=get_quotas)
