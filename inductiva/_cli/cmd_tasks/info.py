"""Prints a task information via CLI."""
from typing import TextIO
import argparse
import sys

from inductiva import tasks, _cli


def task_info(args, fout: TextIO = sys.stdout):
    """Prints a task information."""
    task_id = args.id
    task = tasks.Task(task_id)
    info = task.get_info()
    print(info, file=fout)
    return 0


def register(parser):
    """Register the info tasks command."""
    subparser = parser.add_parser("info",
                                  help="Prints information related to a task.",
                                  formatter_class=argparse.RawTextHelpFormatter)

    subparser.description = ("The `inductiva tasks info` command provides "
                             "an extensive list of information about a task.")

    _cli.utils.add_watch_argument(subparser)
    subparser.add_argument("id",
                           type=str,
                           help="ID of the task to get information about.")
    subparser.set_defaults(func=task_info)
