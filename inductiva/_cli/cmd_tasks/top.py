"""Top command for task."""
from typing import TextIO
import argparse
import sys

from inductiva import tasks


def top(args: argparse.Namespace, fout: TextIO = sys.stdout):
    """Prints the result of the `top -b -H -n 1` command.
    
    This command will list the processes and threads (-H) in batch mode (-b).
    This command will run only once (-n 1) instead of running continuously.
    The result is an instant snapshot of the machine CPU and RAM metrics.
    """
    task_id = args.id
    task = tasks.Task(task_id)
    # pylint: disable=protected-access
    result, return_code = task._top()
    if result:
        print(result, file=fout)
    return return_code


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
