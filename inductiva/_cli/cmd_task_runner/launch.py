"""Launches a Task-Runner via CLI."""
from typing import TextIO
import docker
import argparse
import asyncio
import sys
import os

from inductiva import _cli, constants, _api_key

def generator_thread(queue, generator):
    """Runs a generator in a separate thread."""
    for item in generator:
        queue.put(item)
    queue.put(StopIteration)

def join_generators(*generators):
    """Join multiple generators into one."""
    queue = asyncio.Queue()
    lock = asyncio.Lock()

    for generator in generators:
        yield from generator

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

    for line in join_generators(file_tracker_container.logs(stream=True),
                                task_runner_container.logs(stream=True)):
        print(line.decode("utf-8"), end="", file=fout)

    print(f"File-Tracker launched with container ID: {file_tracker_container.short_id}", file=fout)
    print(f"Task-Runner launched with container ID: {task_runner_container.short_id}", file=fout)


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
