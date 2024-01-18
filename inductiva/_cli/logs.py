"""CLI for logs."""
from inductiva import _cli
from inductiva import logs


def get_task_logs(args):
    """Stream task logs."""
    logs.task_logs.stream_task_logs(args.task_id)


def register_logs_cli(parser):
    _cli.utils.show_help_msg(parser)

    parser.add_argument("task_id",
                        type=int,
                        help="ID of the task to print the logs of.")
    # Register function to call when this subcommand is used
    parser.set_defaults(func=get_task_logs)
