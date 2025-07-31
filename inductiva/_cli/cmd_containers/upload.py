"""
Uploads a Docker image as a converted Apptainer .sif to remote storage.
"""

import argparse
import os
import pathlib
import sys
import tempfile
from inductiva import storage, constants
from inductiva.client.exceptions import ApiException
from inductiva.utils.input_functions import user_confirmation_prompt
from .convert import convert_image
import textwrap


def extract_image_name(image_ref: str) -> str:
    """
    Extracts the base image name from a docker image reference.
    e.g.:
        docker://nginx:alpine → nginx
        nginx:alpine → nginx
        myorg/python → python
        python → python
    """
    # Strip docker:// if present
    image_ref = image_ref.removeprefix("docker://")

    # Split org/image and version tag
    name_part = image_ref.split("/")[-1]
    image_name = "_".join(name_part.split(":"))
    return image_name


def upload_container(args):
    image_name = extract_image_name(args.image)
    default_folder = "my-containers"

    # Handle missing or partial output_path
    output_path = args.output_path

    if not output_path:
        output_path = os.path.join(default_folder, f"{image_name}.sif")
    else:
        output_path = os.path.normpath(output_path)

        # If it's just a filename (no folder), prepend default folder
        if os.sep not in output_path:
            output_path = os.path.join(default_folder, output_path)

        # Ensure .sif extension
        if not output_path.endswith(".sif"):
            output_path += ".sif"

    # Extract folder and filename from the now-final output_path
    folder_name = output_path.split(os.sep)[0]
    filename = os.path.basename(output_path)

    # ---------- pre-flight: does it already exist? ---------------------------
    try:
        contents = storage.listdir(folder_name, print_results=False)
    except ApiException as e:
        if e.status == 404:
            contents = []
        else:
            print(f"Error accessing remote folder '{folder_name}': {e}",
                  file=sys.stderr)
            return

    already_there = any(c["content_name"] == filename for c in contents)
    if already_there:
        if args.overwrite:
            print(f"'{output_path}' exists - will overwrite it "
                  "(because --overwrite).")
            try:
                storage.remove(f"{folder_name}/{filename}")
            except ApiException as e:
                print(f"Failed to delete old file: {e}", file=sys.stderr)
                return
        else:
            confirm = user_confirmation_prompt(
                [filename],
                "File already exists.",
                "File already exists.",
                "Overwrite the existing file?",
                is_all=False,
            )
            if not confirm:
                print("Operation cancelled.")
                return
            storage.remove(f"{folder_name}/{filename}")

    # Create temp folder to hold .sif
    try:
        with tempfile.TemporaryDirectory(dir=constants.TMP_DIR) as tmp_dir:
            sif_folder_path = pathlib.Path(tmp_dir) / folder_name
            os.makedirs(sif_folder_path, exist_ok=True)
            sif_file_path = os.path.join(sif_folder_path, filename)

            convert_args = argparse.Namespace(image=args.image,
                                              output=sif_file_path)

            print(f"Converting {args.image} -> {sif_file_path}...")
            if not convert_image(convert_args):
                print("❌ Conversion failed.")
                return

            print(
                f"Uploading '{sif_folder_path}' to remote dir '{folder_name}'.."
            )
            storage.upload(local_path=sif_folder_path, remote_dir=folder_name)

            # Print the remote path
            print("✅ Upload complete.")
            print("To use the container, instantiate it with:")
            print(f"\t > inductiva://{folder_name}/{filename}")
    except Exception as e:  # pylint: disable=broad-exception-caught
        print("❌ Failed to upload container")
        print(f"Error details: {str(e)}")
        return


def register(parser):
    """Register the upload-container command."""
    subparser = parser.add_parser(
        "upload",
        help="Convert a Docker image to a .sif and upload to remote storage.",
        formatter_class=argparse.RawTextHelpFormatter,
    )

    subparser.description = textwrap.dedent("""\
        Converts a Docker image (from Docker Hub, a local image, or a .tar
        file) into a SIF file using Apptainer, stores it in a temporary folder,
        and uploads that folder to your Inductiva remote storage, making it
        available for use with the Inductiva API.
    """)

    subparser.epilog = textwrap.dedent("""\
        examples:
            # Convert and upload a local Docker image
            $ inductiva containers upload my-simulation-image

            # Convert and upload a Docker Hub CFD image (SU2) with a custom storage path
            $ inductiva containers upload docker://su2code/su2:latest my-containers/su2cfd.sif
    """)

    subparser.add_argument(
        "image",
        type=str,
        help=("Docker image reference. Accepts a:\n"
              "\t- local image name or ID (e.g., python:3.11-slim)\n"
              "\t- Docker Hub reference URL (e.g., docker://nginx:latest)\n"
              "\t- `.tar` archive exported from Docker"),
    )
    subparser.add_argument(
        "output_path",
        nargs="?",
        type=str,
        help=(
            "Optional output path for the `.sif` file in Inductiva remote\n"
            "storage (e.g., `my-containers/nginx.sif`). If omitted, defaults\n"
            "to: `my-containers/<image-name>.sif`."))
    subparser.add_argument(
        "-f",
        "--overwrite",
        action="store_true",
        help="Overwrites the file in remote storage without asking.")

    subparser.set_defaults(func=upload_container)
