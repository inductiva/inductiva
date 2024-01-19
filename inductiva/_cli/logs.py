"""CLI for logs."""
from inductiva import _cli
from inductiva.logs import task_logs


def get_task_logs(args):
    """Stream task logs."""
    streamer = task_logs.TaskLogsStream(args.task_id)
    streamer.stream_task_logs()


def register_logs_cli(parser):
    _cli.utils.show_help_msg(parser)

    parser.add_argument("task_id",
                        type=str,
                        help="ID of the task for which to stream the logs.")
    # Register function to call when this subcommand is used
    parser.set_defaults(func=get_task_logs)
