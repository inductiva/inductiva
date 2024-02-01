"""Kills a tasks by id via CLI."""

import inductiva
from inductiva import constants


def kill_task(args):
    """Kills a task by id."""
    if not args.yes:
        if len(args.id) > constants.MAX_CONFIRMATION_LINES:
            prompt = input(f"You are about to kill {len(args.id)} tasks.\n"
                           "Are you sure you want to proceed (y/[N])? ")
        else:
            print("You are about to kill the following tasks: ")
            for task_id in args.id:
                print(f"  - {task_id}")
            prompt = input("Are you sure you want to proceed (y/[N])? ")

        confirm = prompt.lower() in ["y", "ye", "yes"]
        if not confirm:
            print("Aborted.")
            return

    for task_id in args.id:
        inductiva.tasks.Task(task_id).kill(wait_timeout=args.wait_timeout)


def register(parser):
    """Register the kill task command."""

    subparser = parser.add_parser("kill", help="Kill running tasks.")
    subparser.add_argument("id",
                           type=str,
                           help="ID(s) of the task(s) to kill.",
                           nargs="+")
    subparser.add_argument("-w",
                           "--wait-timeout",
                           type=float,
                           default=None,
                           help="Number of seconds to wait for the kill "
                           "command. If not provided, the system sends "
                           "the request without waiting a response.")
    subparser.add_argument("-y",
                           "--yes",
                           action="store_true",
                           help="Skip kill confirmation.")

    subparser.set_defaults(func=kill_task)
