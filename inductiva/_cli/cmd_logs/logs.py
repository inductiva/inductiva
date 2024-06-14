"""CLI for logs."""
import re
import sys
from typing import Tuple

from .. import utils as cli_utils
from ... import tasks


def _get_task_id_from_mode(mode: str) -> Tuple[bool, str]:
    """Extract the task_id from the mode."""
    match = re.match(r'^(submitted|started)(-(\d+))?$', mode)
    if not match:
        return False, f"Invalid mode format: {mode}"

    status_match = match.group(1)
    offset_match = match.group(3)

    status = status_match
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

    if not cli_utils.is_task_id_valid(mode):
        return False, f"Invalid task_id: {mode}"

    return True, mode


def stream_task_logs(args):
    """Consume the stream logs of a certain task."""
    mode = args.mode.lower()
    result, data = validate_mode_or_task_id(mode)
    if not result:
        print(data, file=sys.stderr)
        return 1

    task_id = data
    task = tasks.Task(data)

    if not task.is_running():
        print(
            f"The current status of task {task_id} is '{task.get_status()}'\n"
            "and the simulation logs are not available for streaming.\n"
            "For more information about the task status, use:\n\n"
            f"  inductiva tasks list --task-id {task_id}\n",
            file=sys.stderr)
        return 1

    consumer = tasks.streams.TaskStreamConsumer(task_id)
    consumer.run_forever()
    return 0


def register(parser):
    cli_utils.show_help_msg(parser)

    parser.add_argument(
        "mode",
        type=str,
        nargs='?',
        default="SUBMITTED",
        help=(
            "Mode of log retrieval. "
            "Use 'SUBMITTED[-n]' for the last submitted tasks, "
            "'STARTED[-n]' for the last started tasks, "
            "or a specific 'task_id' to retrieve logs for a particular task."))

    # Register function to call when this subcommand is used
    parser.set_defaults(func=stream_task_logs)
