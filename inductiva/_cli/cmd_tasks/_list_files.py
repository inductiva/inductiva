"""Command to list directories in a task."""
from typing import TextIO
import argparse
import sys
import asyncio

from inductiva import _cli, tasks


def list_files(args: argparse.Namespace, fout: TextIO = sys.stdout):
    task_id = args.id
    task = tasks.Task(task_id)
    info = task.get_info()
    if info.is_terminal:
        print(
            f"Task {task_id} has terminated.\n"
            "Access its output using:\n\n"
            f"  inductiva tasks download --id {task.id}",
            file=sys.stderr)
        return 1
    if not info.status == "computation-started":
        print(
            f"Task {task_id} has not started yet.\n"
            "Wait for computation to start.",
            file=sys.stderr)
        return 1
    directories = asyncio.run(task._list_files())  # pylint: disable=protected-access
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
