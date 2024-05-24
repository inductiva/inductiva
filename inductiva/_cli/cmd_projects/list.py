"""List the user projetcs information via CLI."""
from collections import defaultdict
from typing import TextIO
import argparse
import sys

import inductiva
from inductiva import _cli
from inductiva.utils import format_utils


def get_projects(_, fout: TextIO = sys.stdout):
    """ Lists the user's projetcs. 

    Lists all the user's projetcs.
    """

    projetcs = inductiva.projects.project.get_projects()
    table = defaultdict(list)
    for p in projetcs:
        table["name"].append(p.name)
        table["nr_tasks"].append(p.num_tasks)

    emph_formatter = format_utils.get_ansi_formatter()

    header_formatters = [
        lambda x: emph_formatter(x.upper(), format_utils.Emphasis.BOLD)
    ]

    table = format_utils.get_tabular_str(table,
                                         header_formatters=header_formatters)

    print(table, file=fout, end="")

    return 0


def register(parser):
    """Register the projetcs list command."""

    subparser = parser.add_parser("list",
                                  help="List the user's projetcs.",
                                  formatter_class=argparse.RawTextHelpFormatter)

    subparser.description = ("The `inductiva projects list` command provides "
                             "an overview of your projects.\n"
                             "It lists all your available projects.")

    _cli.utils.add_watch_argument(subparser)

    subparser.set_defaults(func=get_projects)
