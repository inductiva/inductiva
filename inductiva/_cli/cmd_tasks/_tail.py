"""Tail command for task files."""
from typing import TextIO
import argparse
import sys
import asyncio

from inductiva import _cli, tasks


def tail(args: argparse.Namespace, fout: TextIO = sys.stdout):
    task_id = args.id
    task = tasks.Task(task_id)
    lines = asyncio.run(task._tail_file(args.filename, args.lines))  # pylint: disable=protected-access
    print(lines, file=fout)
    return 0


def register(parser):
    """Register the info tasks command."""
    subparser = parser.add_parser(
        "tail",
        help="Shows the last lines of a file in a task.",
        formatter_class=argparse.RawTextHelpFormatter)

    subparser.description = ("The `inductiva tasks tail` command allows"
                             "to tail a file in a task that is running.")

    _cli.utils.add_watch_argument(subparser)
    subparser.add_argument("id",
                           type=str,
                           help="ID of the task to list directories.")
    subparser.add_argument("filename", type=str, help="File to tail.")
    subparser.add_argument("--lines",
                           "-l",
                           type=int,
                           default=10,
                           help="Number of lines to show.")
    subparser.set_defaults(func=tail)
