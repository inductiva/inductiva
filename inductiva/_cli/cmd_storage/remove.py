"""Remove the user's remote storage contents via CLI."""

import sys

from inductiva import storage


def remove(args):
    """Remove user's remote storage contents."""
    paths = args.path
    confirm = args.confirm
    all_paths = args.all

    if paths and all_paths:
        print("inductiva storage remove: error: "
              "argument path not allowed with argument --all")
        sys.exit(1)

    if not confirm:
        if all_paths:
            prompt = input("Are you sure you want to remove everything "
                           "from your storage(y/[N])")
        else:
            prompt = input(f"Are you sure you want to remove {len(paths)}"
                           " path(s)? (y/[N])")
        confirm = prompt.lower() in ["y", "ye", "yes"]

    if confirm:
        if all_paths:
            storage.rmdir("*", confirm=confirm)
        for path in paths:
            storage.rmdir(path, confirm=confirm)


def register(parser):
    subparser = parser.add_parser("remove",
                                  aliases=["rm"],
                                  help="Remove remote storage entries.")
    subparser.add_argument("path",
                           type=str,
                           nargs="*",
                           help="Path(s) to be removed from remote storage. "
                           "To remove all contents, use \"*\".")
    subparser.add_argument("-y",
                           action="store_true",
                           dest="confirm",
                           default=False,
                           help="Skip confirmation prompt.")
    subparser.add_argument("--all",
                           action="store_true",
                           default=False,
                           help="Remove all contents from remote storage.")

    subparser.set_defaults(func=remove)
