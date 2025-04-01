"""Top command for task."""
import datetime
from typing import TextIO, AsyncGenerator
import argparse
import sys
import asyncio

from inductiva import tasks
from inductiva._cli.cmd_tasks import task_utils
from inductiva.utils import format_utils


def last_modifed_file(args: argparse.Namespace, fout: TextIO = sys.stdout):
    """
    Display the last modified file for a given task.

    This function retrieves and streams information about the most recently 
    modified file associated with a specified task. It validates that the 
    task computation has started before proceeding. If the task is invalid 
    or not started, an error message is printed to `stderr`.
    """
    task_id = args.id
    task = tasks.Task(task_id)
    valid, err_msg = task_utils.validate_task_computation_started(task)
    if not valid:
        print(err_msg, file=sys.stderr)
        return 1
    asyncio.run(stream_task_output(task, fout))
    return 0


async def stream_task_output(task: tasks.Task, fout: TextIO):
    """
    Stream the output of a task's `last_modifed_file` generator to the specified
    output.

    This function gathers and streams the output of the `last_modifed_file`
    method from the given task to the provided file-like object.
    """
    try:
        await asyncio.gather(consume(task.last_modifed_file(), fout))
    except asyncio.CancelledError:
        await task.close_stream()


async def consume(generator: AsyncGenerator, fout: TextIO):
    """
    Consume and write the formatted output from an asynchronous generator to a
    file-like object.

    This function iterates over the provided asynchronous generator, writing 
    each line of output to the specified file-like object.

    Example:
        Most Recent File: /workdir/io9od5da6xh131inmsno0fapm/stdin.txt
        Modification Time: 2025-04-01 09:28:33
        Current Time on Machine: 2025-04-01 09:29:17

        Time Since Last Modification: 0:00:43
    """
    try:
        async for data in generator:

            # Convert timestamps to readable datetime
            most_recent_time = datetime.datetime.fromtimestamp(
                data["most_recent_timestamp"]).strftime("%Y-%m-%d %H:%M:%S")
            now_time = datetime.datetime.fromtimestamp(
                data["now_timestamp"]).strftime("%Y-%m-%d %H:%M:%S")

            # Print the information
            print("", file=fout)
            print(f"Most Recent File: {data["most_recent_file"]}", file=fout)
            print(f"Modification Time: {most_recent_time}", file=fout)
            print(f"Current Time on Machine: {now_time}", file=fout)
            print("", file=fout)
            formatted_seconds = format_utils.seconds_formatter(
                data["time_since_last_mod"])
            print(f"Time Since Last Modification: {formatted_seconds}",
                  file=fout)

    except asyncio.CancelledError:
        pass


def register(parser):
    """Register the last-modifed-file command."""
    subparser = parser.add_parser("last-modified-file",
                                  help="Displays the last modified file.",
                                  formatter_class=argparse.RawTextHelpFormatter)

    subparser.description = (
        "The `inductiva tasks last-modified-file` command displays the last"
        "modified file for the simulation. It also shows the last time the file"
        "was modified and how long it has been since that.")

    subparser.add_argument("id",
                           type=str,
                           help="The ID of the task to get the last modifed "
                           "file.")

    subparser.set_defaults(func=last_modifed_file)
