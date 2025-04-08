"""
List all container files inside a specified (or default) storage folder.
"""

import argparse
from inductiva import storage, client


def list_containers(args):
    default_folder = "my-containers"
    output_path = args.folder

    if not output_path:
        output_path = default_folder

    try:
        storage.listdir(output_path, args.max_results)
    except client.exceptions.ApiException as e:
        if e.status == 404:
            error_msg = f"Folder '{output_path}' not found on remote storage."
            if output_path == default_folder:
                error_msg += ("\nUse the `inductiva containers upload` command "
                              "to upload a container to the remote storage.")
            print(error_msg)
            return False
    except Exception:  # pylint: disable=broad-exception-caught
        print("Unkown error while listing containers.")
        return False


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
        help=(
            "Opt folder path in remote storage to list container files from.\n"
            "If omitted, defaults to 'my-containers'."),
    )

    subparser.add_argument("-m", "--max-results", default=10, type=int)

    subparser.set_defaults(func=list_containers)
