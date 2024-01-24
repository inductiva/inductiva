"""CLI for logs."""
from inductiva import _cli
from inductiva.tasks import streams


def get_task_logs(args):
    """Consume the stream logs of a certain task."""
    consumer = streams.TaskStreamConsumer(args.task_id)
    consumer.run_forever()


def register_logs_cli(parser):
    _cli.utils.show_help_msg(parser)

    parser.add_argument("task_id",
                        type=str,
                        help="ID of the task for which to consume the stream.")
    # Register function to call when this subcommand is used
    parser.set_defaults(func=get_task_logs)
