"""List the user quotas information via CLI."""
from typing import TextIO
import argparse
import sys

from inductiva import users, _cli
from inductiva.utils import format_utils


def get_quotas(_, fout: TextIO = sys.stdout):
    """ Lists the user's quotas. 

    Lists all the user's quotas and the quotas left for the user to use.
    """
    quotas = users.get_quotas()

    name_list = list(quotas)
    max_allowed_list = [value["max_allowed"] for value in quotas.values()]
    in_use_list = [value["in_use"] for value in quotas.values()]

    # Construct the dictionary of lists
    quotas_dict_of_lists = {
        "name": name_list,
        "max_allowed": max_allowed_list,
        "in_use": in_use_list
    }

    emph_formatter = format_utils.get_ansi_formatter()

    header_formatters = [
        lambda x: emph_formatter(x.upper(), format_utils.Emphasis.BOLD)
    ]

    table = format_utils.get_tabular_str(quotas_dict_of_lists,
                                         header_formatters=header_formatters)

    print(table, file=fout, end="")

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
