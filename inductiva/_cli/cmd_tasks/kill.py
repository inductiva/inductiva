"""Kills a tasks by id via CLI."""
import inductiva
from inductiva.utils.input_functions import user_confirmation_prompt


def kill_task(args):
    """Kills a task by id."""
    confirm = args.yes
    if not confirm:
        confirm = user_confirmation_prompt(False, args.id, "",
                                           f"kill {len(args.id)} tasks",
                                           "kill the following tasks")
    if confirm:
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
