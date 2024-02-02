"""Remove the user's remote storage contents via CLI."""

import sys

from inductiva import constants, storage


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
        if all_paths:
            prompt = input("You are about to remove EVERYTHING from your"
                           " remote storage space.\n"
                           "Are you sure you want to proceed (y/[N])? ")
        else:
            if len(paths) > constants.MAX_CONFIRMATION_LINES:
                prompt = input(f"You are about to remove {len(paths)} "
                               "paths from your remote storage space.\n"
                               "Are you sure you want to proceed (y/[N])? ")
            else:
                print("You are about to remove the following paths from "
                      "your remote storage space: ")
                for path in paths:
                    print(f"  - {path}")
                prompt = input("Are you sure you want to proceed (y/[N])? ")
        confirm = prompt.lower() in ["y", "ye", "yes"]

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
