"""Download storage contents via CLI."""

import argparse
import textwrap
from inductiva import storage


def download(args):
    """Download a file or folder from remote storage."""
    storage.download(args.remote_path, args.local_dir, args.decompress)


def register(parser):
    """Register the download command for the user's remote storage."""

    subparser = parser.add_parser(
        "download",
        aliases=["dl"],
        help="Download a file or folder from remote storage.",
        formatter_class=argparse.RawTextHelpFormatter)

    subparser.description = (
        "The `inductiva storage download` command allows you to download a file"
        " or folder from your remote storage, maintaining the original storage "
        "path structure, to your local machine.\n"
        "Specify the remote path, the local destination (optional), and whether"
        " to decompress the content after downloading (default: enabled).\n")
    subparser.add_argument(
        "remote_path",
        type=str,
        help="The path to the file or folder in remote storage to download.",
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
        "--decompress",
        action="store_true",
        help="Decompress the downloaded file or folder if it is compressed.",
    )

    subparser.epilog = textwrap.dedent("""\
        examples:
            # Download and decompress an archived file
            $ inductiva storage download my_data/archive.zip --decompress
    """)

    subparser.set_defaults(func=download)
