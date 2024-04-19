"""List the user quotas information via CLI."""
from typing import TextIO
import argparse
import sys

from inductiva import users, _cli
from inductiva.utils import format_utils


def get_quotas(_, fout: TextIO = sys.stdout):
    """ Lists the user quotas. 

    Lists all the users quotas and the quotas left for the user to use.
    """
    quotas = users.get_quotas()

    # Transpose the list of dictionaries to a dictionary of lists
    quotas = {k: [dic[k] for dic in quotas] for k in quotas[0]}

    emph_formatter = format_utils.get_ansi_formatter()

    header_formatters = [
        lambda x: emph_formatter(x.upper(), format_utils.Emphasis.BOLD)
    ]

    print(format_utils.get_tabular_str(quotas,
                                       headers=quotas.keys(),
                                       header_formatters=header_formatters),
          file=fout,
          end="")

    return 0


def register(parser):
    """Register the user quotas command."""

    subparser = parser.add_parser("quotas",
                                  help="List user quotas.",
                                  formatter_class=argparse.RawTextHelpFormatter)

    subparser.description = ("The `inductiva user quotas` command provides "
                             "an overview of your quotas.\n"
                             "It lists all your available quotas as well as the"
                             "quotas left for you to use\n")

    _cli.utils.add_watch_argument(subparser)

    subparser.set_defaults(func=get_quotas)
