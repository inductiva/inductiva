"""Remove a container file in remote storage."""

import argparse
import sys
from inductiva import storage
from inductiva.utils.input_functions import user_confirmation_prompt


def rm_container(args):
    default_folder = "my-containers"
    folder = args.folder or default_folder
    container_name = args.name

    if not container_name:
        print("Error: Container name is required.", file=sys.stderr)
        return 1

    container_path = f"{folder}/{container_name}"

    # Check if container exists
    try:
        contents = storage.listdir(folder, max_results=100, print_results=False)
    except Exception as e:  # pylint: disable=broad-except
        print(f"Error accessing remote folder '{folder}': {e}", file=sys.stderr)
        return 1

    if not any(
            content["content_name"] == container_name for content in contents):
        print(f"Error: '{container_name}' not found in folder '{folder}'.",
              file=sys.stderr)
        return 1

    # Ask for confirmation if not using --yes
    if not args.confirm:
        confirm = user_confirmation_prompt(
            [container_name],  # items
            "You are about to delete all containers.",
            f"You are about to delete {len([container_name])} containers.",
            "You are about to delete the following container:",
            is_all=False)

        if not confirm:
            print("Operation cancelled.")
            return 0

    try:
        storage.remove_workspace(container_path)
        print(f"✅ Container '{container_name}' removed.")
        return 0
    except Exception as e:  # pylint: disable=broad-except
        print(f"Error removing container: {e}", file=sys.stderr)
        return 1


def register(parser):
    """Register the remove-container command."""
    subparser = parser.add_parser(
        "remove",
        aliases=["rm"],
        help="Remove a container file from remote storage.",
        formatter_class=argparse.RawTextHelpFormatter)

    subparser.description = (
        "Removes a specific container file from a remote storage folder.\n"
        "If no folder is specified, defaults to 'my-containers'.\n"
        "Use the `--yes` flag to skip confirmation prompts.")

    subparser.add_argument(
        "folder",
        nargs="?",
        type=str,
        help=
        "Optional folder path in remote storage. Defaults to 'my-containers'.")

    subparser.add_argument("-n",
                           "--name",
                           required=True,
                           type=str,
                           help="Name of the container to remove.")

    subparser.add_argument(
        "-y",
        "--yes",
        action="store_true",
        dest="confirm",
        default=False,
        help="Skip confirmation prompt and delete immediately.")

    subparser.set_defaults(func=rm_container)
