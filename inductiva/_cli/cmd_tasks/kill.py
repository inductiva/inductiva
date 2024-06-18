"""Kills a tasks by id via CLI."""
import sys

import argparse
import inductiva
from inductiva.tasks.methods import get_all
from inductiva.utils.input_functions import user_confirmation_prompt
from ...localization import translator as __


def kill_task(args):
    """Kills a task by id."""
    kill_all = args.all
    ids = args.id

    if not kill_all and not ids:
        print("No id(s) specified.\n"
              "> Use `inductiva tasks kill -h` for help.")
        return 1

    if ids and kill_all:
        print(
            "inductiva tasks kill: error: "
            "argument id not allowed with argument --all",
            file=sys.stderr)
        return 1
    ids = set(ids)
    if ids:
        tasks = [{"task_id": id} for id in ids]
    else:
        tasks = []
        for status in inductiva.tasks.Task.KILLABLE_STATUSES:
            tasks.extend([{
                "task_id": task.id
            } for task in get_all(status=status)])

    if not tasks:
        print("There are no tasks to kill.")
        return 1

    ids = [task["task_id"] for task in tasks]

    confirm = args.yes or \
        user_confirmation_prompt(ids,
                                __("task-prompt-kill-all"),
                                __("task-prompt-kill-big", len(ids)),
                                __("task-prompt-kill-small"), False
                                )

    if not confirm:
        return 1

    verbosity_level = 1 if kill_all else 2

    for task_id in ids:
        try:
            inductiva.tasks.Task(task_id).kill(wait_timeout=args.wait_timeout,
                                               verbosity_level=verbosity_level)
        except RuntimeError as exc:
            print(f"Error for task {task_id}:", exc)

    return 0


def register(parser):
    """Register the kill task command."""

    subparser = parser.add_parser("kill",
                                  help="Kill running tasks.",
                                  formatter_class=argparse.RawTextHelpFormatter)

    subparser.description = (
        "The `inductiva tasks kill` command terminates specified tasks "
        "on the platform.\n"
        "You can terminate multiple tasks by passive multiple ids.\n"
        "To confirm termination without prompt, use the '-y' or '--yes' "
        "option.\nIf you provide '-w' or '--wait-timeout', the system "
        "does not confirm if the kill command was successful\n")

    subparser.add_argument("id",
                           type=str,
                           help="ID(s) of the task(s) to kill.",
                           nargs="*")
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
    subparser.add_argument("--all",
                           action="store_true",
                           help="Kill all running tasks.")

    subparser.set_defaults(func=kill_task)
