"""Remove the user's remote storage contents via CLI."""

import argparse
from inductiva import storage


def remove(args):
    """Remove user's remote storage contents."""
    path = args.path
    confirm = args.confirm

    if not confirm:
        prompt = input(f"Are you sure you want to remove {path}? (y/n)")
        confirm = prompt.lower() in ["y", "ye", "yes"]

    if confirm:
        print("Removing %s in the remote storage", path)
        storage.rmdir(path, confirm=confirm)


def register(parser):
    subparser = parser.add_parser("remove",
                                  aliases=["rm"],
                                  help="Remove remote storage entries.",
                                  formatter_class=argparse.RawTextHelpFormatter)
    subparser.description = (
        "The `inductiva storage remove` command deletes specified data from the platform.\n"
        "It targets a specific path for removal. Use with caution as this action is irreversible.\n\n"
        "Use `inductiva storage remove \"*\"` to remove all folders in your storage."
    )
    subparser.add_argument("path",
                           type=str,
                           help="Path to be removed from remote storage. "
                           "To remove all contents, use \"*\".")
    subparser.add_argument("-y",
                           "--yes",
                           action="store_true",
                           dest="confirm",
                           help="Skip remove confirmation.",
                           default=False)
    subparser.set_defaults(func=remove)
