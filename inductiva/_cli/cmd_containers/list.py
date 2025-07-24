"""
List all container files inside a specified (or default) storage folder.
"""

import argparse
from inductiva import storage, client
import textwrap


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
        help=("List all container files in remote storage, "
              "including their size and estimated cost."),
        formatter_class=argparse.RawTextHelpFormatter)

    subparser.description = textwrap.dedent("""\
        The `inductiva containers list` command lists all container files in
        remote storage under the default containers folder (`my-containers/`),
        including their size and estimated cost.
        
        You can also specify other container folders to list their contents.
    """)

    subparser.add_argument(
        "folder",
        nargs="?",
        type=str,
        help="Path to a containers folder in remote storage.",
    )

    subparser.add_argument(
        "--max-results",
        "-m",
        default=10,
        type=int,
        help="Maximum number of container files to list.",
    )

    subparser.epilog = textwrap.dedent("""\
        examples:
            $ inductiva containers list
            NAME             SIZE        CREATION TIME
            container1.sif   200.00 MB   26/03, 16:41:14
            container2.sif   100.00 MB   26/03, 16:41:14
            container3.sif   300.00 MB   26/03, 16:41:14

            Total storage size used:
                    Volume: 0.59 GB
                    Cost: 0.02 US$/month
    """)

    subparser.set_defaults(func=list_containers)
