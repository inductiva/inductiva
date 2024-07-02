"""CLI for logs."""
import json
import re
import sys
from typing import Tuple

from .. import utils as cli_utils
from ... import tasks, ApiException


def _get_task_id_from_mode(mode: str) -> Tuple[bool, str]:
    """Extract the task_id from the mode."""
    match = re.match(r"^(submitted|started)(-(\d+))?$", mode)
    if not match:
        return False, f"Invalid mode format: {mode}"

    offset_match = match.group(3)

    status = match.group(1)
    offset = (int(offset_match) + 1) if offset_match else 1

    # If the status is submitted, send None to the API to get all the tasks.
    # The endpoint already returns the tasks in the submitted order by default
    if status == "submitted":
        status = None

    task_list = tasks.get(last_n=offset, status=status)

    if not task_list:
        return False, f"No '{status}' tasks found."

    return True, task_list[-1].id


def validate_mode_or_task_id(mode: str) -> Tuple[bool, str]:
    """Validate the mode or task_id."""
    if mode.startswith("submitted") or mode.startswith("started"):
        return _get_task_id_from_mode(mode)

    return True, mode


def _check_if_task_is_running(task: tasks.Task) -> Tuple[bool, str]:
    """Check if the task is running.
    
    This is done inside a try-except block to catch any exceptions that may
    occur when trying to get the task status. Mainly if the task does not exist.
    """

    is_running = False
    try:
        info = task.get_info()
        is_running = info.is_running
    except ApiException as e:
        return False, f"{ json.loads(e.body)['detail']}"

    if not is_running:
        return False, (
            f"The current status of task {task.id} is '{info.status}'\n"
            "and the simulation logs are not available for streaming.\n"
            "For more information about the task status, use:\n\n"
            f"  inductiva tasks list --id {task.id}\n")

    return True, task.id


def stream_task_logs(args):
    """Consume the stream logs of a certain task."""
    mode = args.mode.lower()
    result, data = validate_mode_or_task_id(mode)
    if not result:
        print(data, file=sys.stderr)
        return 1

    io_stream = None if (
        args.stdout and args.stderr
    ) else "std_out" if args.stdout else "std_err" if args.stderr else None

    no_color = True if io_stream else args.no_color

    task_id = data
    task = tasks.Task(data)

    result, data = _check_if_task_is_running(task)
    if not result:
        print(data, file=sys.stderr)
        return 1

    consumer = tasks.streams.TaskStreamConsumer(task_id,
                                                io_stream=io_stream,
                                                no_color=no_color)
    consumer.run_forever()
    return 0


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
    # Register function to call when this subcommand is used
    parser.set_defaults(func=stream_task_logs)
