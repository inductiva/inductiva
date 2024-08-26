"""Remove the user's remote storage contents via CLI."""
import argparse

import inductiva
from inductiva.utils.input_functions import user_confirmation_prompt
from ...localization import translator as __


def remove(args):
    """Remove user's remote storage contents."""
    paths = args.id
    confirm = args.confirm

    paths = set(paths)
    if not confirm:
        confirm = user_confirmation_prompt(
            paths, __("storage-prompt-remove-all"),
            __("storage-prompt-remove-big", len(paths)),
            __("storage-prompt-remove-small"), False)

    if confirm:
        for task_id in paths:
            inductiva.tasks.Task(task_id).remove_remote_files()

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
                           nargs="+",
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

    subparser.set_defaults(func=remove)
