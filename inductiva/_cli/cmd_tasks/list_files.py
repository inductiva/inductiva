"""Command to list directories in a task."""
from typing import TextIO
import argparse
import sys

from inductiva import _cli, tasks


def list_files(args: argparse.Namespace, fout: TextIO = sys.stdout):
    task_id = args.id
    task = tasks.Task(task_id)

    directories, ret_code = task.list_files()
    if directories:
        print(directories, file=fout)
    return ret_code


def register(parser):
    """Register the info tasks command."""
    subparser = parser.add_parser(
        "list-files",
        help="Lists the current files of a running task.",
        formatter_class=argparse.RawTextHelpFormatter)

    subparser.description = ("The `inductiva tasks list-files` command lists "
                             "the contents of a task's working directory while "
                             "the task is in progress. (Experimental)")

    _cli.utils.add_watch_argument(subparser)
    subparser.add_argument("id", type=str, help="ID of the task to list files.")
    subparser.set_defaults(func=list_files)
