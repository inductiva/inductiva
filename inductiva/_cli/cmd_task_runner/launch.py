"""Launches a Task-Runner via CLI."""
from typing import TextIO
import argparse
import threading
import sys
import os
import platform
from inductiva import _cli, constants, _api_key, api_url
from inductiva.resources import machine_groups

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


def launch_task_runner_gcp(args):
    """Launches a Task-Runner on GCP.
    
    Note: BYOC (Bring Your Own Cloud) machine groups can only have a single VM.
    """
    machine_groups.MachineGroup(
        machine_type=args.machine_type,
        zone=args.zone,
        provider="GCP",
        spot=args.spot,
        mg_name=args.machine_group_name,
        byoc=True,
        max_idle_time=args.max_idle_time,
    ).start()


def launch_task_runner(args, fout: TextIO = sys.stdout):
    """Launches a Task-Runner."""
    if args.provider == "gcp":
        launch_task_runner_gcp(args)
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
        "  gcp    - Launch on Google Cloud Platform")

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

    gcp_group.add_argument("--spot",
                           action="store_true",
                           help="Use spot instance.")

    gcp_group.add_argument("--max-idle-time",
                           type=int,
                           help="Idle time in minutes before "
                                "auto-termination (default: 3).")

    subparser.set_defaults(func=launch_task_runner)
