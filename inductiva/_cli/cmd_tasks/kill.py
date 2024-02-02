"""Kills a tasks by id via CLI."""
import inductiva
from inductiva.utils.input_functions import user_confirmation_prompt
from ...localization import translator as __


def kill_task(args):
    """Kills a task by id."""
    confirm = args.yes or user_confirmation_prompt(
        args.id, __("user-prompt-kill-all"),
        __("user-prompt-kill-big", len(args.id)), __("user-prompt-kill-small"),
        False)

    if confirm:
        for task_id in args.id:
            inductiva.tasks.Task(task_id).kill(wait_timeout=args.wait_timeout)
    return 0


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
