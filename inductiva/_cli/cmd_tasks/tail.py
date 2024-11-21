from typing import TextIO
import argparse
import sys
import asyncio
import logging

logging.getLogger().setLevel(level=logging.WARNING)

from inductiva import _cli, tasks


def tail(args: argparse.Namespace, fout: TextIO = sys.stdout):
    task_id = args.id
    task = tasks.Task(task_id)
    lines = asyncio.run(
        task.file_operation(tasks.Operations.TAIL, filename=args.filename))
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
    subparser.set_defaults(func=tail)
