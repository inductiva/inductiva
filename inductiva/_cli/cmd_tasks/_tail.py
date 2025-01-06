"""Tail command for task files."""
from typing import TextIO
import argparse
import sys
import asyncio

from inductiva import _cli, tasks


def tail(args: argparse.Namespace, fout: TextIO = sys.stdout):
    task_id = args.id
    task = tasks.Task(task_id)
    asyncio.run(consume(task, args, fout))
    return 0


async def consume(task: tasks.Task, args: argparse.Namespace, fout: TextIO):
    prefix = "\033[s"
    try:
        async for lines in task._tail_file(  # pylint: disable=protected-access
                args.filename, args.lines, args.follow):
            print(f"{prefix}{lines}", file=fout, end="", flush=True)
            prefix = "\033[u\033[s"
    except asyncio.CancelledError:
        await task.close_stream()


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
    subparser.add_argument(
        "--follow",
        "-f",
        action="store_true",
        help="Keep the file open and show new lines.",
    )
    subparser.set_defaults(func=tail)
