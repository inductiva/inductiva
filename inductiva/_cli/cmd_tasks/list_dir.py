from typing import TextIO
import argparse
import sys
import asyncio
import logging

logging.getLogger().setLevel(level=logging.WARNING)

from inductiva import _cli, tasks


def list_directories(args: argparse.Namespace, fout: TextIO = sys.stdout):
    task_id = args.id
    task = tasks.Task(task_id)
    directories = asyncio.run(task.file_operation(tasks.Operations.LIST))
    print(directories, file=fout)
    return 0


def register(parser):
    """Register the info tasks command."""
    subparser = parser.add_parser(
        "ls",
        help="Lists directories in a task recursively.",
        formatter_class=argparse.RawTextHelpFormatter)

    subparser.description = ("The `inductiva tasks ls` command provides "
                             "the directories generated during a task.")

    _cli.utils.add_watch_argument(subparser)
    subparser.add_argument("id",
                           type=str,
                           help="ID of the task to list directories.")
    subparser.set_defaults(func=list_directories)
