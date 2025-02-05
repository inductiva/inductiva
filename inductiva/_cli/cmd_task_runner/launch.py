"""Launches a Task-Runner via CLI."""
from typing import TextIO
import argparse
import threading
import sys
import os
import platform
from inductiva import _cli, constants, _api_key, api_url

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


def launch_task_runner(args, fout: TextIO = sys.stdout):
    """Launches a Task-Runner."""
    if not _docker_imported:
        print(
            "Docker Python API not installed, please run "
            "'pip install inductiva[task-runner]' to install it",
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

    task_runner_container = client.containers.run(
        image=constants.TASK_RUNNER_IMAGE,
        name="task-runner",
        environment={
            "USER_API_KEY": _api_key.get(),
            "API_URL": api_url,
            "MACHINE_GROUP_NAME": args.machine_group_name,
            "HOST_NAME": args.hostname,
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
        "a way to launch a Task-Runner on the platform.")

    _cli.utils.add_watch_argument(subparser)
    subparser.add_argument(
        "machine_group_name",
        type=str,
        help="Name of the machine group to launch the Task-Runner.")

    subparser.add_argument("--hostname",
                           "-ho",
                           type=str,
                           default=platform.uname().node,
                           help="Hostname of the Task-Runner.")

    subparser.add_argument("--detach",
                           "-d",
                           action="store_true",
                           help="Run the task-runner in the background.")

    subparser.set_defaults(func=launch_task_runner)
