"""Top command for task."""
from typing import TextIO, AsyncGenerator
import argparse
import sys
import asyncio

from inductiva import tasks
from inductiva._cli.cmd_tasks import task_utils


def top(args: argparse.Namespace, fout: TextIO = sys.stdout):
    """Prints the result of the `top -b -H -n 1` command.
    
    This command will list the processes and threads (-H) in batch mode (-b).
    This command will run only once (-n 1) instead of running continuously.
    The result is an instant snapshot of the machine CPU and RAM metrics."
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
    Stream the output of a task's `list_top` generator to the specified output.

    This function gathers and streams the output of the `list_top` method 
    from the given task to the provided file-like object.
    """
    try:
        await asyncio.gather(consume(task.list_top(), fout))
    except asyncio.CancelledError:
        await task.close_stream()


async def consume(generator: AsyncGenerator, fout: TextIO):
    """
    Consume and write the output from an asynchronous generator to a file-like
    object.

    This function iterates over the provided asynchronous generator, writing 
    each line of output to the specified file-like object.
    """
    try:
        async for lines in generator:
            print(lines, file=fout, end="", flush=True)
    except asyncio.CancelledError:
        pass


def register(parser):
    """Register the `top` command for monitoring task resource usage."""
    subparser = parser.add_parser(
        "top",
        help="Displays the output of the `top` command from the machine running"
        " the task.",
        formatter_class=argparse.RawTextHelpFormatter)

    subparser.description = (
        "The `inductiva tasks top` command streams the output of the"
        "`top -b -H -n 1` command executed on the machine where the "
        "task is running. This allows real-time monitoring of system "
        "processes and resource usage for the task.")

    subparser.add_argument("id",
                           type=str,
                           help="The ID of the task to monitor.")

    subparser.set_defaults(func=top)
