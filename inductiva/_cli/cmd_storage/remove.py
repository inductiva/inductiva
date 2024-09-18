"""Remove the user's remote storage contents via CLI."""
import argparse
import sys

import inductiva
from inductiva.tasks.methods import get_all
from ...localization import translator as __
from inductiva.utils.input_functions import user_confirmation_prompt


def remove(args):
    """Remove user's remote storage contents."""
    is_all = args.all
    task_ids = args.id
    confirm = args.confirm

    if is_all and task_ids:
        print(
            "inductiva storage remove: error: "
            "argument id(s) not allowed with argument --all/-a",
            file=sys.stderr)
        return 1
    if not is_all and not task_ids:
        print(
            "inductiva storage remove: error: "
            "argument id(s) or --all/-a required",
            file=sys.stderr)
        return 1

    task_ids = set(task_ids)

    if is_all:
        all_tasks = get_all()
        task_ids = [task.id for task in all_tasks]

    if not confirm:
        confirm = user_confirmation_prompt(
            task_ids, __("storage-prompt-remove-all"),
            __("storage-prompt-remove-big", len(task_ids)),
            __("storage-prompt-remove-small"), is_all)

    if not confirm:
        return 0

    failed = False
    for task_id in task_ids:
        if not inductiva.tasks.Task(task_id).remove_remote_files(verbose=False):
            print(f"Failed to remove the following task storage: {task_id}\n",
                  file=sys.stderr)
            failed = True

    if failed:
        return 1

    print("All tasks storage removed successfully.")
    return 0


def register(parser):
    subparser = parser.add_parser("remove",
                                  aliases=["rm"],
                                  help="Remove remote storage entries.",
                                  formatter_class=argparse.RawTextHelpFormatter)
    subparser.description = (
        "The `inductiva storage remove` command deletes specified data"
        " from the platform.\n"
        "It targets a specific task files for removal. Use with caution as "
        "this action is irreversible.\n\n")
    subparser.add_argument("id",
                           type=str,
                           nargs="*",
                           help="Id(s) of the task(s) to remove from remote "
                           "storage.")

    subparser.add_argument(
        "-y",
        "--yes",
        action="store_true",
        dest="confirm",
        default=False,
        help="Sets any confirmation values to \"yes\" "
        "automatically. Users will not be asked for "
        "confirmation to remove path(s) from remote storage.")

    subparser.add_argument("-a",
                           "--all",
                           action="store_true",
                           default=False,
                           help="Remove all tasks from remote storage.")

    subparser.set_defaults(func=remove)
