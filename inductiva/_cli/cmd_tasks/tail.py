"""Tail command for task files."""
from typing import TextIO, AsyncGenerator
import argparse
import sys
import asyncio

from inductiva import tasks
from inductiva._cli.cmd_tasks import task_utils


def tail(args: argparse.Namespace, fout: TextIO = sys.stdout):
    task_id = args.id
    task = tasks.Task(task_id)
    valid, err_msg = task_utils.validate_task_computation_started(task)
    if not valid:
        print(err_msg, file=sys.stderr)
        return 1
    asyncio.run(gather_tasks(task, args, fout))
    return 0


async def gather_tasks(task: tasks.Task, args: argparse.Namespace,
                       fout: TextIO):
    generators = [
        task.tail_file(filename, args.lines, args.follow)
        for filename in args.filename
    ]
    tail_tasks = [
        asyncio.create_task(consume(generator, fout))
        for generator in generators
    ]
    try:
        await asyncio.gather(*tail_tasks)
    except asyncio.CancelledError:
        for tail_task in tail_tasks:
            tail_task.cancel()
        await task.close_stream()


async def consume(generator: AsyncGenerator, fout: TextIO):
    try:
        async for lines in generator:
            print(lines, file=fout, end="", flush=True)
    except asyncio.CancelledError:
        pass


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
