"""Tail command for task files."""
from typing import TextIO
import argparse
import sys
import asyncio

from inductiva import tasks


def tail(args: argparse.Namespace, fout: TextIO = sys.stdout):
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
    asyncio.run(consume(task, args, fout))
    return 0


async def consume(task: tasks.Task, args: argparse.Namespace, fout: TextIO):
    try:
        async for lines in task._tail_file(  # pylint: disable=protected-access
                args.filename, args.lines, args.follow):
            print(lines, file=fout, end="", flush=True)
    except asyncio.CancelledError:
        await task.close_stream()


def register(parser):
    """Register the info tasks command."""
    subparser = parser.add_parser("tail",
                                  formatter_class=argparse.RawTextHelpFormatter)

    subparser.description = (
        "The `inductiva tasks tail` shows the last lines "
        "of a file in the directory of a task that is running. "
        "(Experimental)")

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
