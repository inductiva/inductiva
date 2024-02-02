"""Remove the user's remote storage contents via CLI."""

import sys

from inductiva import storage
from inductiva.utils.input_functions import user_confirmation_prompt
from ...localization import translator as __


def remove(args):
    """Remove user's remote storage contents."""
    paths = args.path
    confirm = args.confirm
    all_paths = args.all

    if paths and all_paths:
        print(
            "inductiva storage remove: error: "
            "argument path not allowed with argument --all",
            file=sys.stderr)
        sys.exit(1)

    if not confirm:
        confirm = user_confirmation_prompt(
            paths, __("user-prompt-remove-all"),
            __("user-prompt-remove-big", len(paths)),
            __("user-prompt-remove-small"), all_paths)

    if confirm:
        if all_paths:
            storage.rmdir("*", confirm=True)
        for path in paths:
            storage.rmdir(path, confirm=True)


def register(parser):
    subparser = parser.add_parser("remove",
                                  aliases=["rm"],
                                  help="Remove remote storage entries.")
    subparser.add_argument("path",
                           type=str,
                           nargs="*",
                           help="Path(s) to be removed from remote storage. "
                           "To remove all contents, use \"*\".")
    subparser.add_argument(
        "-y",
        action="store_true",
        dest="confirm",
        default=False,
        help="Sets any confirmation values to \"yes\" "
        "automatically. Users will not be asked for "
        "confirmation to remove path(s) from remote storage.")
    subparser.add_argument("--all",
                           action="store_true",
                           default=False,
                           help="Remove all contents from remote storage.")

    subparser.set_defaults(func=remove)
