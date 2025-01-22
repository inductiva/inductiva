"""Removes the Task-Runner via CLI."""
from typing import TextIO
import argparse
import sys
from inductiva import _cli

_docker_imported = True
try:
    import docker
    from docker.errors import DockerException
except ImportError:
    _docker_imported = False


def remove_task_runner(_, fout: TextIO = sys.stdout):
    """Removes the Task-Runner."""
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

    if not (task_runner or file_tracker):
        print("Task-Runner is already stopped", file=fout)
        return

    if file_tracker:
        print("Removing the File-Tracker container.", file=fout)
        file_tracker[0].stop()
    if task_runner:
        print("Removing the Task-Runner container.", file=fout)
        task_runner[0].stop()


def register(parser):
    """Register the remove task-runner command."""
    subparser = parser.add_parser("remove",
                                  help="Removes the Task-Runner.",
                                  formatter_class=argparse.RawTextHelpFormatter)

    subparser.description = (
        "The `inductiva task-runner remove` command provides "
        "a way to remove the Task-Runner running in the background.")

    _cli.utils.add_watch_argument(subparser)
    subparser.set_defaults(func=remove_task_runner)
