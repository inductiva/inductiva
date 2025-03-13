"""Download storage contents via CLI."""

import argparse
from inductiva import storage


def download(args):
    """Download a file or folder from remote storage."""
    storage.download(args.remote_path, args.local_dir, args.uncompress)


def register(parser):
    """Register the download command for the user's remote storage."""

    subparser = parser.add_parser(
        "download",
        aliases=["dl"],
        help="Download a file or folder from remote storage.",
        formatter_class=argparse.RawTextHelpFormatter)

    subparser.description = (
        "The `inductiva storage download` command allows you to download files "
        "or folders from your remote storage.\n"
        "Specify the remote path, the local destination (optional), and whether"
        " to uncompress the content after downloading (default: enabled).\n")
    subparser.add_argument(
        "remote_path",
        type=str,
        help="The remote path to the file or folder to download.",
    )
    subparser.add_argument(
        "local_dir",
        type=str,
        nargs="?",
        default="",
        help=("The local directory where the downloaded content will be saved. "
              "Defaults to the current working directory if not specified."),
    )
    subparser.add_argument(
        "-u",
        "--uncompress",
        action="store_true",
        help=("Uncompress the downloaded file or folder if it is compressed. "
              "This option is enabled by default."),
    )

    subparser.set_defaults(func=download)
