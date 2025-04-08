"""
Uploads a Docker image as a converted Apptainer .sif to remote storage.
"""

import argparse
import os
import tempfile
from inductiva import storage
from .convert import convert_image


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
    image_name = name_part.split(":")[0]
    return image_name


def upload_container(args):
    image_name = extract_image_name(args.image)
    default_folder = "my-containers"

    # Handle missing or partial output_path
    output_path = args.output_path

    if not output_path:
        output_path = f"{default_folder}/{image_name}.sif"
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

    # Create temp folder to hold .sif
    try:
        with tempfile.TemporaryDirectory() as tmp_dir:
            sif_folder_path = os.path.join(tmp_dir, folder_name)
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
    subparser.description = (
        "Converts a Docker image (from Docker Hub or local) into a .sif file "
        "using Apptainer,\n stores it temporarily in a folder, and uploads that"
        " folder to the system's remote storage.")

    subparser.add_argument(
        "image",
        type=str,
        help="Docker image reference (e.g., python:3.11-slim or docker://...).",
    )
    subparser.add_argument(
        "output_path",
        nargs="?",
        type=str,
        help=(
            "Optional output path for the .sif file, my-containers/nginx.sif.\n"
            "If omitted, defaults to my-containers/<image-name>.sif."),
    )

    subparser.set_defaults(func=upload_container)
