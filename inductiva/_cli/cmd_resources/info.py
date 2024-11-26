"""Prints a machine group information via CLI."""
from typing import TextIO
import argparse
import sys

from inductiva import resources, _cli

def machine_group_info(args, fout: TextIO = sys.stdout):
    """Prints a task information."""
    machine_group_name = args.name
    machine_group = resources.machine_groups.get_by_name(machine_group_name)
    return 0


def register(parser):
    """Register the info machine group command."""
    subparser = parser.add_parser("info",
                                  help="Prints information related to a machine group.",
                                  formatter_class=argparse.RawTextHelpFormatter)

    subparser.description = ("The `inductiva resources info` command provides "
                             "information about a machine group.")

    _cli.utils.add_watch_argument(subparser)
    subparser.add_argument("name",
                           type=str,
                           help="name of the machine group to get information about.")
    subparser.set_defaults(func=machine_group_info)
