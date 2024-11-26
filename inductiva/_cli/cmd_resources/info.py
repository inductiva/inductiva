"""Prints a machine group information via CLI."""
from typing import TextIO
import argparse
import sys

from inductiva import resources, _cli
from inductiva.utils import format_utils


def machine_group_info(args, fout: TextIO = sys.stdout):
    """Prints a task information."""
    machine_group_name = args.name
    machine_group = resources.machine_groups.get_by_name(machine_group_name)
    machines = machine_group.__dict__["machines"]

    print(f"Showing machines of machine group: {machine_group_name}", file=fout)

    table_dict = resources.machine_groups.get_machine_dict(machines)

    formatters = {
        "Started": [format_utils.datetime_formatter,],
        "Last seen": [format_utils.datetime_formatter,],
    }

    emph_formatter = format_utils.get_ansi_formatter()

    header_formatters = [
        lambda x: emph_formatter(x.upper(), format_utils.Emphasis.BOLD)
    ]

    print(format_utils.get_tabular_str(table_dict,
                                       formatters=formatters,
                                       header_formatters=header_formatters),
          file=fout,
          end="")

    return 0


def register(parser):
    """Register the info machine group command."""
    subparser = parser.add_parser(
        "info",
        help="Prints information related to a machine group.",
        formatter_class=argparse.RawTextHelpFormatter)

    subparser.description = ("The `inductiva resources info` command provides "
                             "information about a machine group.")

    _cli.utils.add_watch_argument(subparser)
    subparser.add_argument(
        "name",
        type=str,
        help="name of the machine group to get information about.")
    subparser.set_defaults(func=machine_group_info)
