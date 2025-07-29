"""Kills a tasks by id via CLI."""
import sys

import argparse
import textwrap
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
        # pylint: disable=protected-access
        for status in inductiva.tasks.Task._KILLABLE_STATUSES:
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

    subparser.description = textwrap.dedent("""\
        The `inductiva tasks kill` command terminates the specified tasks
        on the platform. You can terminate multiple tasks by passing multiple
        task IDs. To confirm termination without a prompt, use the `-y`
        or `--yes` option. If you provide `-w` or `--wait-timeout`, the system
        does not confirm whether the termination was successful.

        Note: The `inductiva tasks kill` command does not stop the machine
        where the task is running. It only terminates the task itself, leaving
        the computational resources active and available to run other tasks.
    """)

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

    subparser.epilog = textwrap.dedent("""\
        examples:
            $ inductiva tasks kill cmvsc9qhz5iy86f6pef8uyxqt
            You are about to kill the following tasks:
                - cmvsc9qhz5iy86f6pef8uyxqt 
            Are you sure you want to proceed (y/[N])? y
            Successfully sent kill request for task cmvsc9qhz5iy86f6pef8uyxqt.
    """)

    subparser.set_defaults(func=kill_task)
