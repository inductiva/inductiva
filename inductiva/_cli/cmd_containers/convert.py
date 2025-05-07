"""
Converts a Docker image to an Apptainer .sif file via CLI.
"""

import argparse
import threading
import sys
import os
import tempfile
from typing import TextIO
import uuid

from inductiva import constants

try:
    import docker
    from docker.errors import DockerException, ImageNotFound, NotFound
except ImportError:
    docker = None


def generator_thread(container, fout: TextIO):
    """Stream logs from a container."""
    for line in container.logs(stream=True):
        print(line.decode("utf-8"), end="", file=fout)


def join_container_streams(*containers, fout: TextIO = sys.stdout):
    """Join multiple container log streams into one output."""
    container_threads = []
    for container in containers:
        thread = threading.Thread(target=generator_thread,
                                  args=(container, fout))
        thread.start()
        container_threads.append(thread)
    for thread in container_threads:
        thread.join()


def convert_image(args, fout: TextIO = sys.stdout):
    """
    Converts a Docker image to an Apptainer .sif file.
    
    The source image may be provided as:
      - a Docker Hub URL (e.g. docker://python:3.11-slim),
      - a local Docker image reference (by image name or ID),
      - or a path to a tar archive.
    
    The output file will be saved to the specified location.
    """
    if docker is None:
        print(
            "Docker Python API not installed. Please run "
            "'pip install inductiva[task-runner]' to install it.",
            file=fout,
        )
        return False

    try:
        client = docker.from_env()
    except DockerException as e:
        print(f"Failed to connect to Docker: {e}", file=fout)
        return False

    # Prepare output file path
    output = os.path.abspath(args.output)
    outdir = os.path.dirname(output)
    if outdir and not os.path.exists(outdir):
        os.makedirs(outdir, exist_ok=True)

    source = args.image
    conversion_source = None
    tmp_tar_path = None

    # Determine the conversion source
    if source.startswith("docker://"):
        # Use the docker:// URL directly.
        conversion_source = source
    elif os.path.exists(source):
        # Assume it's a tar file.
        conversion_source = f"docker-archive:///{os.path.abspath(source)}"
    else:
        # Assume it's a local Docker image reference.
        try:
            image_obj = client.images.get(source)
        except ImageNotFound:
            print(
                f"Local image '{source}' not found. "
                "\nThe image must be pulled or built in the system",
                file=fout)
            return False

        # Save the image to a temporary tar file.
        with tempfile.NamedTemporaryFile(suffix=".tar",
                                         delete=False,
                                         dir=constants.TMP_DIR) as tmp:
            tmp_tar_path = tmp.name
            with open(tmp_tar_path, "wb") as f:
                for chunk in image_obj.save(named=True):
                    f.write(chunk)
        conversion_source = f"docker-archive:///{tmp_tar_path}"

    mounts = []

    # If .tar or local docker image
    if tmp_tar_path:
        tar_host_path = os.path.abspath(tmp_tar_path)
        tar_host_dir = os.path.dirname(tar_host_path)
        tar_filename = os.path.basename(tar_host_path)
        tar_inside_container = f"/input/{tar_filename}"
        conversion_source = f"docker-archive://{tar_inside_container}"
        mounts.append(
            docker.types.Mount(target="/input",
                               source=tar_host_dir,
                               type="bind"))

    output_filename = os.path.basename(output)
    out_mount_path = "/output"
    mounts.append(
        docker.types.Mount(target=out_mount_path,
                           source=os.path.abspath(outdir),
                           type="bind"))

    command = [
        "build", f"{out_mount_path}/{output_filename}", conversion_source
    ]

    container = None
    try:
        container = client.containers.run(
            image=constants.TASK_RUNNER_IMAGE,
            name=f"apptainer-converter-{uuid.uuid4().hex[:8]}",
            entrypoint="/usr/bin/apptainer",
            command=command,
            mounts=mounts,
            platform="linux/amd64",
            detach=True,
            auto_remove=False,
        )

        print("Conversion started with container ID:",
              container.short_id,
              file=fout)
        join_container_streams(container, fout=fout)
        print(f"Conversion complete. Output saved at: {output}", file=fout)

    except KeyboardInterrupt:
        print("\n⚠️ Conversion interrupted by user. Cleaning up...", file=fout)
        if container:
            try:
                container.kill()
                container.remove(force=True)
                print("🧹 Container stopped and removed.", file=fout)
            except Exception as cleanup_err:  # pylint: disable=broad-exception-caught
                print(f"⚠️ Failed to clean up container: {cleanup_err}",
                      file=fout)
        return False

    except Exception as e:  # pylint: disable=broad-exception-caught
        print(f"Error during conversion: {e}", file=fout)
        return False

    finally:
        # Clean up the container if it was created.
        try:
            if container:
                container.remove(force=True)
        except NotFound:
            pass

    # Clean up the temporary tar file if one was created.
    if tmp_tar_path and os.path.exists(tmp_tar_path):
        try:
            os.remove(tmp_tar_path)
        except Exception as e:  # pylint: disable=broad-exception-caught
            print(
                f"Warning: could not remove temporary file {tmp_tar_path}: {e}",
                file=fout)
            return False

    return True


def register(parser):
    """Register the convert command."""
    subparser = parser.add_parser(
        "convert",
        help="Convert a Docker image to an Apptainer .sif file.",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    subparser.description = (
        "The `inductiva containers convert` command converts a Docker "
        "image—specified as a Docker Hub URL,\n a local Docker image reference,"
        " or a .tar archive—into an Apptainer .sif file using\n the Apptainer.")
    subparser.add_argument(
        "image",
        type=str,
        help=("Docker image reference. A docker:// URL, a local image name/ID, "
              "or a path to a .tar file."),
    )
    subparser.add_argument(
        "output",
        type=str,
        help="Path to save the resulting .sif file.",
    )
    subparser.set_defaults(func=convert_image)
