"""Tail command for task files."""
from typing import TextIO
import argparse
import sys

from inductiva import tasks


def tail(args: argparse.Namespace, fout: TextIO = sys.stdout):
    task_id = args.id
    task = tasks.Task(task_id)

    return task.tail_files(args.filename, args.lines, args.follow, fout)


def register(parser):
    """Register the info tasks command."""
    subparser = parser.add_parser(
        "tail",
        help="Prints the last -l/--lines of a file of a running task.",
        formatter_class=argparse.RawTextHelpFormatter)

    subparser.description = (
        "The `inductiva tasks tail` shows the last lines "
        "of a file in the directory of a task that is running. "
        "(Experimental)")

    subparser.add_argument("id",
                           type=str,
                           help="ID of the task to list directories.")
    subparser.add_argument("filename",
                           type=str,
                           nargs="+",
                           help="File to tail.")
    subparser.add_argument("--lines",
                           "-l",
                           type=int,
                           default=10,
                           help="Number of lines to show.")
    subparser.add_argument(
        "--follow",
        "-f",
        action="store_true",
        help="Keep the file open and show new lines.",
    )
    subparser.set_defaults(func=tail)
