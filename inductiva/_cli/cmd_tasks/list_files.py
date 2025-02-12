"""Command to list directories in a task."""
from typing import TextIO
import argparse
import sys
import asyncio

from inductiva import _cli, tasks
from inductiva._cli.cmd_tasks import task_utils


def list_files(args: argparse.Namespace, fout: TextIO = sys.stdout):
    task_id = args.id
    task = tasks.Task(task_id)
    valid, err_msg = task_utils.validate_task_computation_started(task)
    if not valid:
        print(err_msg, file=sys.stderr)
        return 1
    directories = asyncio.run(task.list_files())  # pylint: disable=protected-access
    print(directories, file=fout)
    return 0


def register(parser):
    """Register the info tasks command."""
    subparser = parser.add_parser("list-files",
                                  formatter_class=argparse.RawTextHelpFormatter)

    subparser.description = ("The `inductiva tasks list-files` command lists "
                             "the contents of a task's working directory while "
                             "the task is in progress. (Experimental)")

    _cli.utils.add_watch_argument(subparser)
    subparser.add_argument("id", type=str, help="ID of the task to list files.")
    subparser.set_defaults(func=list_files)
