"""
List all container files inside a specified (or default) storage folder.
"""

import argparse
from inductiva import storage


def list_containers(args):
    default_folder = "my-containers"
    output_path = args.folder

    if not output_path:
        output_path = default_folder

    storage.listdir(output_path, args.max_results)


def register(parser):
    """Register the upload-container command."""
    subparser = parser.add_parser(
        "list",
        aliases=["ls"],
        help="List all container files in remote storage.",
        formatter_class=argparse.RawTextHelpFormatter)

    subparser.description = (
        "Lists all container files in the folder within remote storage.\n"
        "If no folder is provided, defaults to the 'my-containers' directory.")

    subparser.add_argument(
        "folder",
        nargs="?",
        type=str,
        help=
        ("Optional folder path in remote storage to list container files from.\n"
         "If omitted, defaults to 'my-containers'."),
    )

    subparser.add_argument("-m", "--max-results", default=10, type=int)

    subparser.set_defaults(func=list_containers)
