"""Kills a tasks by id via CLI."""

import inductiva


def kill_task(args):
    """Kills a task by id."""
    inductiva.tasks.Task(args.id).kill(wait_timeout=args.wait_timeout)


def register(parser):
    """Register the kill task command."""

    subparser = parser.add_parser("kill", help="Kill a running task.")
    subparser.add_argument("id", type=str, help="ID(s) of the task to kill.", nargs='+')
    subparser.add_argument("-w",
                           "--wait-timeout",
                           type=float,
                           default=None,
                           help="Timeout to wait for the kill to complete.")

    subparser.set_defaults(func=kill_task)
