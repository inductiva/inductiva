"""Launches a Task-Runner via CLI."""
from typing import TextIO
import argparse
import threading
import sys
import os
import platform
import subprocess
import tempfile
from inductiva import _cli, constants, _api_key, api_url
from inductiva.resources import byoc_gcp

_docker_imported = True
try:
    import docker
    from docker.errors import DockerException
except ImportError:
    _docker_imported = False


def generator_thread(container, fout: TextIO):
    """Runs a generator in a separate thread."""
    for line in container.logs(stream=True):
        print(line.decode("utf-8"), end="", file=fout)


def join_container_streams(*containers, fout: TextIO = sys.stdout):
    """Join multiple containers into one."""
    container_threads = []

    for container in containers:
        thread = threading.Thread(target=generator_thread,
                                  args=(container, fout))
        thread.start()
        container_threads.append(thread)

    for thread in container_threads:
        thread.join()


def check_gcloud_installed() -> bool:
    """Check if gcloud CLI is installed and authenticated."""
    try:
        subprocess.run(["gcloud", "--version"],
                       capture_output=True,
                       text=True,
                       check=True)
        return True
    except (subprocess.CalledProcessError, FileNotFoundError):
        return False


def check_gcloud_auth() -> bool:
    """Check if gcloud is authenticated."""
    try:
        result = subprocess.run(
            ["gcloud", "auth", "list", "--filter=status:ACTIVE"],
            capture_output=True,
            text=True,
            check=True)
        return "ACTIVE" in result.stdout
    except subprocess.CalledProcessError:
        return False


def launch_task_runner_gcp(args, fout: TextIO = sys.stdout):
    """Launches a Task-Runner on GCP."""

    if not check_gcloud_installed():
        print("Error: gcloud CLI is not installed or not in PATH.", file=fout)
        print(
            "Please install gcloud CLI: "
            "https://cloud.google.com/sdk/docs/install",
            file=fout)
        print("Or install with: pip install 'inductiva[gcp]'", file=fout)
        return

    if not check_gcloud_auth():
        print("Error: gcloud is not authenticated.", file=fout)
        print("Please run 'gcloud auth login' to authenticate.", file=fout)
        return

    api_key = _api_key.get()
    if not api_key:
        print("Error: No API key found. Please set your API key first.",
              file=fout)
        return

    startup_script = byoc_gcp.create_gcp_startup_script()

    with tempfile.NamedTemporaryFile(mode="w", suffix=".sh", delete=False) as f:
        f.write(startup_script)
        script_path = f.name

    try:
        cmd = [
            "gcloud", "compute", "instances", "create", args.machine_group_name,
            "--zone", args.zone, "--machine-type", args.machine_type,
            "--image-family", args.image_family, "--image-project",
            args.image_project, "--scopes",
            "https://www.googleapis.com/auth/cloud-platform", "--metadata",
            f"INDUCTIVA_API_KEY={api_key},"
            f"INDUCTIVA_API_URL={api_url},"
            f"MACHINE_GROUP_NAME={args.machine_group_name}",
            "--metadata-from-file", f"startup-script={script_path}"
        ]

        if args.preemptible:
            cmd.append("--preemptible")

        if args.hostname:
            cmd.extend(["--metadata",
                       f"TASK_RUNNER_HOSTNAME={args.hostname}"])

        print(
            f"Creating GCP VM '{args.machine_group_name}' "
            f"in zone '{args.zone}'...",
            file=fout)
        print(f"Machine type: {args.machine_type}", file=fout)
        if args.preemptible:
            print("Using preemptible instance (spot pricing)", file=fout)

        result = subprocess.run(cmd, capture_output=True, text=True,
                               check=False)

        if result.returncode == 0:
            print("GCP VM created successfully!", file=fout)
            print(f"VM Name: {args.machine_group_name}", file=fout)
            print(f"Zone: {args.zone}", file=fout)
            print(
                "The task-runner will start automatically "
                "once the VM is ready.",
                file=fout)
        else:
            print("Failed to create GCP VM:", file=fout)
            print(result.stderr, file=fout)

    except (subprocess.CalledProcessError, OSError) as e:
        print(f"Error creating GCP VM: {e}", file=fout)
    finally:
        try:
            os.unlink(script_path)
        except OSError:
            pass


def launch_task_runner(args, fout: TextIO = sys.stdout):
    """Launches a Task-Runner."""
    if args.provider == "gcp":
        launch_task_runner_gcp(args, fout)
    else:  # local
        launch_task_runner_local(args, fout)


def launch_task_runner_local(args, fout: TextIO = sys.stdout):
    """Launches a Task-Runner locally using Docker."""
    if not _docker_imported:
        print(
            "Docker Python API not installed, please run "
            "'pip install \"inductiva[task-runner]\"' to install it",
            file=fout)
        return

    try:
        client = docker.from_env()
    except DockerException as e:
        print(f"Failed to connect to Docker: {e}", file=fout)
        print(
            "Please make sure Docker is running "
            "and you have the necessary permissions.",
            file=fout)
        return

    task_runner = client.containers.list(filters={"name": "task-runner"})
    file_tracker = client.containers.list(filters={"name": "file-tracker"})

    if task_runner or file_tracker:
        print(
            "Task-Runner already running. "
            "Please stop it before launching a new one.",
            file=fout)
        return

    client.images.pull(constants.FILE_TRACKER_IMAGE, platform="linux/amd64")
    file_tracker_container = client.containers.run(
        image=constants.FILE_TRACKER_IMAGE,
        name="file-tracker",
        environment={
            "API_URL": api_url,
            "USER_API_KEY": _api_key.get(),
        },
        volumes={
            "workdir": {
                "bind": "/workdir",
                "mode": "rw"
            },
        },
        network="host",
        platform="linux/amd64",
        detach=True,
        auto_remove=True,
    )

    client.images.pull(constants.TASK_RUNNER_IMAGE, platform="linux/amd64")

    apptainer_path = "apptainer"
    os.makedirs(apptainer_path, exist_ok=True)
    os.chmod(apptainer_path, 0o777)
    apptainer_full_path = os.path.abspath(apptainer_path)

    hostname = args.hostname or platform.uname().node

    task_runner_container = client.containers.run(
        image=constants.TASK_RUNNER_IMAGE,
        name="task-runner",
        environment={
            "USER_API_KEY": _api_key.get(),
            "API_URL": api_url,
            "MACHINE_GROUP_NAME": args.machine_group_name,
            "HOST_NAME": hostname,
        },
        mounts=[
            docker.types.Mount(target="/executer-images",
                               source=apptainer_full_path,
                               type="bind")
        ],
        volumes={
            "workdir": {
                "bind": "/workdir",
                "mode": "rw"
            },
        },
        network="host",
        privileged=True,
        platform="linux/amd64",
        detach=True,
        auto_remove=True,
    )

    print(
        "File-Tracker launched "
        f"with container ID: {file_tracker_container.short_id}",
        file=fout)
    print(
        "Task-Runner launched "
        f"with container ID: {task_runner_container.short_id}",
        file=fout)

    if not args.detach:
        try:
            join_container_streams(task_runner_container,
                                   file_tracker_container)
        except KeyboardInterrupt:
            print("Interrupted. Stopping containers...", file=fout)
            file_tracker_container.stop()
            task_runner_container.stop()


def register(parser):
    """Register the launch task-runner command."""
    subparser = parser.add_parser("launch",
                                  help="Launches a Task-Runner.",
                                  formatter_class=argparse.RawTextHelpFormatter)

    subparser.description = (
        "The `inductiva task-runner launch` command provides "
        "a way to launch a Task-Runner on different providers.\n\n"
        "Providers:\n"
        "  local  - Launch locally using Docker (default)\n"
        "  gcp    - Launch on Google Cloud Platform\n\n"
        "Local Provider:\n"
        "  - Pulls required Docker images\n"
        "  - Creates necessary volumes and directories\n"
        "  - Launches file-tracker and task-runner containers\n"
        "  - Monitors containers (unless --detach is used)\n\n"
        "GCP Provider:\n"
        "  - Creates a GCP compute instance\n"
        "  - Installs Docker and pulls required images\n"
        "  - Launches file-tracker and task-runner containers\n"
        "  - Monitors task-runner and auto-deletes VM when it stops\n\n"
        "Prerequisites:\n"
        "  Local: Docker installed and running\n"
        "  GCP:   gcloud CLI installed and authenticated\n"
        "  Both:  Valid Inductiva API key configured")

    _cli.utils.add_watch_argument(subparser)

    subparser.add_argument(
        "--provider",
        "-p",
        choices=["local", "gcp"],
        default="local",
        help="Provider to use for launching the task-runner (default: local).")

    subparser.add_argument(
        "machine_group_name",
        type=str,
        help="Name of the machine group to launch the Task-Runner.")

    subparser.add_argument("--hostname",
                           "-ho",
                           type=str,
                           help="Hostname of the Task-Runner.")

    subparser.add_argument(
        "--detach",
        "-d",
        action="store_true",
        help="Run the task-runner in the background (local only).")

    # GCP-specific arguments
    gcp_group = subparser.add_argument_group("GCP-specific options")
    gcp_group.add_argument(
        "--zone",
        "-z",
        type=str,
        default="europe-west1-b",
        help="GCP zone where the VM will be created (default: europe-west1-b).")

    gcp_group.add_argument("--machine-type",
                           "-t",
                           type=str,
                           default="c2d-standard-8",
                           help="GCP machine type (default: c2d-standard-8).")

    gcp_group.add_argument("--image-family",
                           type=str,
                           default="ubuntu-2204-lts",
                           help="GCP image family (default: ubuntu-2204-lts).")

    gcp_group.add_argument("--image-project",
                           type=str,
                           default="ubuntu-os-cloud",
                           help="GCP image project (default: ubuntu-os-cloud).")

    gcp_group.add_argument(
        "--preemptible",
        action="store_true",
        help="Use preemptible instance (spot pricing) for cost savings.")

    subparser.set_defaults(func=launch_task_runner)
