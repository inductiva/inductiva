"""Remove the user's remote storage contents via CLI."""

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
                                  description="Remove remote storage entries.",
                                  help="Remove remote storage entries.")
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
