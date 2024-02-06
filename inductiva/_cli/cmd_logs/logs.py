"""CLI for logs."""
import sys

from .. import utils as cli_utils
from ... import tasks


def stream_task_logs(args):
    """Consume the stream logs of a certain task."""
    task_id = args.task_id
    task = tasks.Task(task_id)

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

    parser.add_argument("task_id",
                        type=str,
                        help="ID of the task for which to consume the stream.")

    # Register function to call when this subcommand is used
    parser.set_defaults(func=stream_task_logs)
