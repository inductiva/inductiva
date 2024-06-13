"""CLI for logs."""
import sys

from .. import utils as cli_utils
from ... import tasks


def stream_task_logs(args):
    """Consume the stream logs of a certain task."""

    if args.stdout and args.stderr:
        print(
            "inductiva logs: error: "
            "argument --stdout not allowed with argument --stderr",
            file=sys.stderr)
        return 1

    io_stream = "std_out" if args.stdout else "std_err" if args.stderr else None
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

    consumer = tasks.streams.TaskStreamConsumer(task_id, io_stream=io_stream)
    consumer.run_forever()
    return 0


def register(parser):
    cli_utils.show_help_msg(parser)

    parser.add_argument("task_id",
                        type=str,
                        help="ID of the task for which to consume the stream.")
    parser.add_argument("--stdout",
                        action="store_true",
                        help="Consumes the standard output stream of the task.")
    parser.add_argument("--stderr",
                        action="store_true",
                        help="Consumes the standard error stream of the task.")

    # Register function to call when this subcommand is used
    parser.set_defaults(func=stream_task_logs)
