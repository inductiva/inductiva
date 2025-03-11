"""Download storage contents via CLI."""

import argparse
from inductiva import storage


def download(args):
    """Download a file or folder from remote storage."""
    storage.download(args.remote_path, args.local_dir, args.uncompress)


def register(parser):
    """Register the download command for the user's remote storage."""

    subparser = parser.add_parser("download",
                                  aliases=["dl"],
                                  help="Download storage file or folder.",
                                  formatter_class=argparse.RawTextHelpFormatter)

    subparser.description = (
        "The `inductiva storage download` command allows you to download files "
        "or folders from your remote storage.\n"
        "Specify the remote path, the local destination, and whether or not to "
        "uncompress the content after downloading.\n")
    subparser.add_argument("remote_path",
                           type=str,
                           help="The remote path to the file or folder.")
    subparser.add_argument("local_dir",
                           type=str,
                           default="",
                           help="Local directory to save the content.")
    subparser.add_argument("-u",
                           "--uncompress",
                           action="store_true",
                           help="Uncompress the content after download.")

    subparser.set_defaults(func=download)
