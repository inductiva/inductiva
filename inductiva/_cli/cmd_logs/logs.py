"""Consume the stream of logs of a task."""

from inductiva.tasks import streams


def stream_task_logs(args):
    """Consume the stream logs of a certain task."""
    consumer = streams.TaskStreamConsumer(args.task_id)
    consumer.run_forever()


def register(parser):
    """Register the logs commands to stream the logs of a user's tasks."""
    parser.add_argument("task_id",
                        type=str,
                        help="ID of the task for which to consume the stream.")

    parser.set_defaults(func=stream_task_logs)
