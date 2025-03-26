"""
Uploads a Docker image as a converted Apptainer .sif to remote storage.
"""

import argparse
import os
import tempfile
from inductiva import storage
from inductiva import _cli


def upload_container(args):
    output_subpath = args.output_path
    folder_name = os.path.normpath(output_subpath).split(os.sep)[0]
    filename = os.path.basename(output_subpath)

    # 2. Create temp dir to hold folder and .sif file
    with tempfile.TemporaryDirectory() as tmp_dir:
        sif_folder_path = os.path.join(tmp_dir, folder_name)
        os.makedirs(sif_folder_path, exist_ok=True)
        sif_file_path = os.path.join(sif_folder_path, filename)

        convert_args = argparse.Namespace(image=args.image,
                                          output=sif_file_path)

        print(f"Converting {args.image} -> {sif_file_path}...")
        _cli.convert_image(convert_args)

        # 4. Upload entire folder to storage
        print(f"Uploading '{sif_folder_path}' to remote dir '{folder_name}'...")
        storage.upload(local_path=sif_folder_path, remote_dir=folder_name)

        print("✅ Upload complete.")


def register(parser):
    """Register the upload-container command."""
    subparser = parser.add_parser(
        "upload-container",
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
        type=str,
        help=("Path (folder/filename) to save and upload the .sif file, "
              "e.g., my-containers/python.sif."),
    )

    subparser.set_defaults(func=upload_container)
