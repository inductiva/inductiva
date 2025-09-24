"""Remove the user's remote storage contents via CLI."""
import argparse
import textwrap

import inductiva
from ...localization import translator as __
from inductiva.utils.input_functions import user_confirmation_prompt


def remove(args):
    """Remove user's remote storage contents."""
    if not args.confirm:
        confirm = user_confirmation_prompt(
            args.paths,
            __("storage-prompt-remove-all"),
            __("storage-prompt-remove-big", len(args.paths)),
            __("storage-prompt-remove-small"),
            args.all,
        )
        if not confirm:
            return

    if args.all:
        root_dir = inductiva.storage.listdir(max_results=None)
        args.paths = [directory["content_name"] for directory in root_dir]

    for path in args.paths:
        inductiva.storage.remove(remote_path=path)


def register(parser):
    subparser = parser.add_parser(
        "remove",
        aliases=["rm"],
        help="Remove files or directories from remote storage.",
        formatter_class=argparse.RawTextHelpFormatter)
    subparser.description = (
        "The `inductiva storage remove` command deletes specific files or"
        "folders from your Inductiva remote storage.\n"
        "Use it with caution. This action is irreversible and will permanently "
        "remove the selected files or directories from your remote storage.\n")

    subparser.add_argument("paths",
                           type=str,
                           nargs="*",
                           help="Remote path(s) to remove from storage.")
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
                           help="Remove all data from remote storage.")

    subparser.epilog = textwrap.dedent("""\
        examples:
            $ inductiva storage remove 0bet8jrpp2gz974n42nsd9n2p/
            You are about to remove the following paths from your remote storage space:
            - 0bet8jrpp2gz974n42nsd9n2p/
            Are you sure you want to proceed (y/[N])? y
            Removing '0bet8jrpp2gz974n42nsd9n2p/' from remote storage...
            Successfully removed '0bet8jrpp2gz974n42nsd9n2p/' from remote storage.
    """)

    subparser.set_defaults(func=remove)
