"""Launches a Task-Runner via CLI."""
from typing import TextIO
import argparse
import threading
import sys
import os
import platform
import subprocess
import signal
from inductiva import _cli, constants, _api_key, api_url

_docker_imported = True
try:
    import udocker
except ImportError:
    _docker_imported = False


def generator_thread(container, fout: TextIO):
    """Runs a generator in a separate thread."""
    for line in container.stdout:
        print(line, file=fout)


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

    # try:
    #     client = docker.from_env()
    # except DockerException as e:
    #     print(f"Failed to connect to Docker: {e}", file=fout)
    #     print(
    #         "Please make sure Docker is running "
    #         "and you have the necessary permissions.",
    #         file=fout)
    #     return

    # task_runner = client.containers.list(filters={"name": "task-runner"})
    # file_tracker = client.containers.list(filters={"name": "file-tracker"})

    # if task_runner or file_tracker:
    #     print(
    #         "Task-Runner already running. "
    #         "Please stop it before launching a new one.",
    #         file=fout)
    #     return
    workdir_path = "workdir"
    os.makedirs(workdir_path, exist_ok=True)
    workdir_abs_path = os.path.abspath(workdir_path)

    subprocess.run(["udocker", "pull", "--platform", "linux/amd64", constants.FILE_TRACKER_IMAGE])
    subprocess.run(["udocker", "setup", "--execmode=S1", constants.FILE_TRACKER_IMAGE])

    cmd = f"""
    udocker run \
    --env='API_URL={api_url}' \
    --env='USER_API_KEY={_api_key.get()}' \
    --name=file-tracker \
    --volume={workdir_abs_path}:/workdir \
    --platform=linux/amd64 \
    --rm {constants.FILE_TRACKER_IMAGE}
    """

    file_tracker = subprocess.Popen(cmd,
                                    stdout=subprocess.PIPE, shell=True)
    
    # client.images.pull(constants.FILE_TRACKER_IMAGE, platform="linux/amd64")
    # file_tracker_container = client.containers.run(
    #     image=constants.FILE_TRACKER_IMAGE,
    #     name="file-tracker",
    #     environment={
    #         "API_URL": api_url,
    #         "USER_API_KEY": _api_key.get(),
    #     },
    #     volumes={
    #         "workdir": {
    #             "bind": "/workdir",
    #             "mode": "rw"
    #         },
    #     },
    #     network="host",
    #     platform="linux/amd64",
    #     detach=True,
    #     auto_remove=True,
    # )
    apptainer_path = "apptainer"
    os.makedirs(apptainer_path, exist_ok=True)
    os.chmod(apptainer_path, 0o777)
    apptainer_full_path = os.path.abspath(apptainer_path)

    subprocess.run(["udocker", "pull", "--platform", "linux/amd64",  constants.TASK_RUNNER_IMAGE])
    subprocess.run(["udocker", "setup", "--execmode=S1", constants.TASK_RUNNER_IMAGE])

    cmd = f"""
    udocker run \
    --env='API_URL={api_url}' \
    --env='USER_API_KEY={_api_key.get()}' \
    --env='MACHINE_GROUP_NAME={args.machine_group_name}' \
    --env='HOST_NAME={args.hostname}' \
    --name=task-runner \
    --volume={workdir_abs_path}:/workdir \
    --volume={apptainer_full_path}:/executer-images \
    --volume=/etc/hosts:/etc/hosts \
    --platform=linux/amd64 \
    --rm {constants.TASK_RUNNER_IMAGE}
    """
                                     
    task_runner = subprocess.Popen(cmd,
                                    stdout=subprocess.PIPE, shell=True)

    # client.images.pull(constants.TASK_RUNNER_IMAGE, platform="linux/amd64")

    # task_runner_container = client.containers.run(
    #     image=constants.TASK_RUNNER_IMAGE,
    #     name="task-runner",
    #     environment={
    #         "USER_API_KEY": _api_key.get(),
    #         "API_URL": api_url,
    #         "MACHINE_GROUP_NAME": args.machine_group_name,
    #         "HOST_NAME": args.hostname,
    #     },
    #     mounts=[
    #         docker.types.Mount(target="/executer-images",
    #                            source=apptainer_full_path,
    #                            type="bind")
    #     ],
    #     volumes={
    #         "workdir": {
    #             "bind": "/workdir",
    #             "mode": "rw"
    #         },
    #     },
    #     network="host",
    #     privileged=True,
    #     platform="linux/amd64",
    #     detach=True,
    #     auto_remove=True,
    # )

    # print(
    #     "File-Tracker launched "
    #     f"with container ID: {file_tracker_container.short_id}",
    #     file=fout)
    # print(
    #     "Task-Runner launched "
    #     f"with container ID: {task_runner_container.short_id}",
    #     file=fout)

    # if not args.detach:
    #     try:
    #         join_container_streams(task_runner_container,
    #                                file_tracker_container)
    #     except KeyboardInterrupt:
    #         print("Interrupted. Stopping containers...", file=fout)
    #         file_tracker_container.stop()
    #         task_runner_container.stop()

    try:
        join_container_streams(task_runner,
                                file_tracker)
    except KeyboardInterrupt:
        print("Interrupted. Stopping containers...", file=fout)
        os.kill(file_tracker.pid, signal.SIGKILL)
        os.kill(task_runner.pid, signal.SIGKILL)
        subprocess.run(["udocker", "rm", "task-runner"])
        subprocess.run(["udocker", "rm", "file-tracker"])

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
