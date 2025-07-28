"""Remove a container file in remote storage."""

import argparse
import sys
from inductiva import storage
from inductiva.utils.input_functions import user_confirmation_prompt
import textwrap


def rm_container(args):
    folder = args.folder
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
        storage.remove(container_path)
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
        help="Remove a container file from your Inductiva remote storage.",
        formatter_class=argparse.RawTextHelpFormatter)

    subparser.description = textwrap.dedent("""\
        The `inductiva containers remove` command removes a specific container
        file from your Inductiva remote storage. If no folder is specified,
        defaults to `my-containers/`.
    
        Use the flag `--yes` to skip confirmation prompts.
        
        This action is irreversible and should be used with caution.
    """)

    subparser.add_argument("folder",
                           nargs="?",
                           type=str,
                           default="my-containers/",
                           help="Path to folder in remote storage.")

    subparser.add_argument("-n",
                           "--name",
                           required=True,
                           type=str,
                           help="Name of the container file to remove.")

    subparser.add_argument(
        "-y",
        "--yes",
        action="store_true",
        dest="confirm",
        default=False,
        help="Skip confirmation prompt and delete immediately.")

    subparser.epilog = textwrap.dedent("""\
        examples:
            # Remove a container with confirmation
            $ inductiva containers rm -n nginx.sif

            # Remove a container without confirmation prompt
            $ inductiva containers rm -n nginx.sif -y

            # Remove a container from a specific folder
            $ inductiva containers rm my-custom-folder -n my-container.sif
    """)

    subparser.set_defaults(func=rm_container)
