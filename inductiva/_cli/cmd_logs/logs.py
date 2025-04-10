"""CLI for logs."""
import sys

from .. import utils as cli_utils
from ... import tasks


def stream_task_logs_tail(args):
    task_id = args.mode.lower()
    task = tasks.Task(task_id)

    filename = "stdout.txt" if args.stdout else "stderr.txt"
    files = ["stdout.txt", "stderr.txt"
            ] if args.stdout == args.stderr else [filename]

    return task.tail_files(files, 10, True, sys.stdout)


def register(parser):
    cli_utils.show_help_msg(parser)

    parser.add_argument(
        "mode",
        type=str,
        nargs="?",
        default="SUBMITTED",
        help=
        ("Mode of log retrieval. "
         "Use 'SUBMITTED' for the last submitted task, "
         "'SUBMITTED-1' for the second last submitted task, and so on. "
         "'STARTED' and 'STARTED-n' follow the same pattern for started tasks. "
         "Or, use a specific 'task_id' to retrieve logs for a particular task."
        ))

    parser.add_argument("--stdout",
                        action="store_true",
                        help="Consumes the standard output stream of the task.")
    parser.add_argument("--stderr",
                        action="store_true",
                        help="Consumes the standard error stream of the task.")
    parser.add_argument("--no-color",
                        action="store_true",
                        help="Disables the colorized output.")
    parser.add_argument(
        "--wait",
        "-w",
        action="store_true",
        help=("Wait for the task to start running before consuming the logs.\n"
              "Without this flag, the logs will be consumed immediately\nor "
              "returns an error if the task is not running."))
    # Register function to call when this subcommand is used
    parser.set_defaults(func=stream_task_logs_tail)
