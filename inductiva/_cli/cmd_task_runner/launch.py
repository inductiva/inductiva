"""Launches a Task-Runner via CLI."""
from typing import TextIO
import docker
import argparse
import asyncio
import queue
import threading
import sys
import os

from inductiva import _cli, constants, _api_key

def generator_thread(container, fout: TextIO):
    """Runs a generator in a separate thread."""
    for line in container.logs(stream=True):
        print(line.decode("utf-8"), end="", file=fout)


def join_container_streams(*containers, fout: TextIO = sys.stdout):
    """Join multiple containers into one."""
    container_threads = []

    for container in containers:
        thread = threading.Thread(target=generator_thread, args=(container, fout))
        thread.start()
        container_threads.append(thread)
    
    print("Terminating threads...")
    for thread in container_threads:
        thread.join()

def launch_task_runner(args, fout: TextIO = sys.stdout):
    """Launches a Task-Runner."""
    client = docker.from_env()
    file_tracker_container = client.containers.run(
        image=constants.FILE_TRACKER_IMAGE,
        environment={
            "USER_API_KEY": _api_key.get(),
        },
        volumes={
            "workdir": {'bind': '/workdir', 'mode': 'rw'},
        },
        network="host",
        detach=True,
    )
    task_runner_container = client.containers.run(
        image=constants.TASK_RUNNER_IMAGE,
        environment={
            "USER_API_KEY": _api_key.get(),
            "MACHINE_GROUP_NAME": args.machine_group_name,
            "HOST_NAME": args.hostname,
        },
        volumes={
            "workdir": {'bind': '/workdir', 'mode': 'rw'},
            "apptainer": {'bind': '/executer-images', 'mode': 'rw'},
        },
        network="host",
        detach=True,
    )

    print(f"File-Tracker launched with container ID: {file_tracker_container.short_id}", file=fout)
    print(f"Task-Runner launched with container ID: {task_runner_container.short_id}", file=fout)

    try:
        join_container_streams(task_runner_container, file_tracker_container)
    except KeyboardInterrupt:
        print("Interrupted. Stopping containers...", file=fout)
        task_runner_container.stop()
        file_tracker_container.stop()

def register(parser):
    """Register the launch task-runner command."""
    subparser = parser.add_parser("launch",
                                  help="Launches a Task-Runner.",
                                  formatter_class=argparse.RawTextHelpFormatter)

    subparser.description = ("The `inductiva task-runner launch` command provides "
                             "a way to launch a Task-Runner on the platform.")

    _cli.utils.add_watch_argument(subparser)
    subparser.add_argument("machine_group_name",
                           type=str,
                           help="Name of the machine group to launch the Task-Runner.")
    
    subparser.add_argument("--hostname",
                            "-ho",
                            type=str,
                            default=os.uname().nodename,
                            help="Hostname of the Task-Runner.")
    

    subparser.set_defaults(func=launch_task_runner)
