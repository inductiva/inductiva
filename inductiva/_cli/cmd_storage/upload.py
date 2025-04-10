"""Upload storage contents via CLI."""

import argparse
from inductiva import storage


def upload(args):
    """Upload a file or folder to remote storage."""
    storage.upload(args.local_path, args.remote_dir)


def register(parser):
    """Register the upload command for the user's remote storage."""

    subparser = parser.add_parser(
        "upload",
        aliases=["ul"],
        help="Upload a file or folder to remote storage.",
        formatter_class=argparse.RawTextHelpFormatter)

    subparser.description = (
        "The `inductiva storage upload` command allows you to upload files "
        "or folders to your remote storage.\n"
        "Specify the local path and the remote destination.\n\n"
        "Example:\n"
        "    inductiva storage upload local/path/file_or_directory remote_dir\n"
    )

    subparser.add_argument(
        "local_path",
        type=str,
        help="The local path to the file or folder to upload.",
    )
    subparser.add_argument(
        "remote_dir",
        type=str,
        help=(
            "The remote directory where the uploaded content will be stored."),
    )

    subparser.set_defaults(func=upload)
