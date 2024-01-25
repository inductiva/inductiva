"""Manage running/completed tasks on the Inductiva API."""
import pathlib
import contextlib
import time
import json
from absl import logging
from typing import Dict, Any, List, Optional
from typing_extensions import TypedDict
import datetime
from dateutil import parser

import inductiva
from inductiva.client import models
from inductiva import api, types
from inductiva.client.apis.tags import tasks_api
from inductiva.utils import files, format_utils, data, output_contents

import warnings

_TASK_TERMINAL_STATUSES = {
    models.TaskStatusCode.SUCCESS, models.TaskStatusCode.FAILED,
    models.TaskStatusCode.KILLED, models.TaskStatusCode.EXECUTERFAILED,
    models.TaskStatusCode.EXECUTERTERMINATED,
    models.TaskStatusCode.SPOTINSTANCEPREEMPTED, models.TaskStatusCode.ZOMBIE
}


class Task:
    """Represents a running/completed task on the Inductiva API.

    Example usage:
        task = scenario.simulate(...)
        final_status = task.wait()
        info = task.get_info() # dictionary with info about the task
        task.download_outputs(
            filenames=["file1.txt", "file2.dat"] # download only these files
        )

    Attributes:
        id: The task ID.
        _api: Instance of TasksApi (from the generated API client).
    """

    def __init__(self, task_id: str):
        """Initialize the instance from a task ID."""
        self.id = task_id
        self._api = tasks_api.TasksApi(api.get_client())
        self._info = None
        self._status = None

    @classmethod
    def from_api_info(cls, info: Dict[str, Any]) -> "Task":

        task = cls(info["task_id"])
        task._info = info
        task._status = models.TaskStatusCode(info["status"])

        # TODO(luispcunha): construct correct output class from API info.

        return task

    @contextlib.contextmanager
    def sync_context(self):
        """Enter context manager for blocking sync execution.

        This turns an asynchronous task into a blocking task.
        If an exception/ctrl+c is caught while in the context manager, the
        remote task is killed.

        Usage:
            task = scenario.simulate(...)

            with task.sync_context():
                # If an exception happens here or ctrl+c is pressed, the task
                # will be killed.
                task.wait()

        """
        try:
            yield
        except KeyboardInterrupt:
            logging.info("Caught SIGINT: terminating blocking task...")
            self.kill()
        except Exception as e:
            logging.error("Caught exception: terminating blocking task...")
            self.kill()
            raise e

    def get_status(self) -> models.TaskStatusCode:
        """Get status of the task.

        This method issues a request to the API.
        """
        # If the task is in a terminal status and we already have the status,
        # return it without refreshing it from the API.
        if self._status is not None and self._status in _TASK_TERMINAL_STATUSES:
            return self._status

        resp = self._api.get_task_status(self._get_path_params())

        return models.TaskStatusCode(resp.body["status"])

    def get_info(self) -> Dict[str, Any]:
        """Get a dictionary with information about the task.

        Information includes e.g., "task_id", "status", timestamps
        ("create_time", "input_submit_time, "start_time", "end_time"),
        among others.

        This method issues a request to the API.
        """
        # If the task is in a terminal status and we already have the info,
        # return it without refreshing it from the API.
        if self._info is not None and self._status in _TASK_TERMINAL_STATUSES:
            return self._info

        params = self._get_path_params()
        resp = self._api.get_task(params, skip_deserialization=True).response

        info = json.loads(resp.data.decode("utf-8"))
        status = models.TaskStatusCode(info["status"])

        self._info = info
        self._status = status

        return info

    def wait(self, polling_period: int = 5) -> models.TaskStatusCode:
        """Wait for the task to complete.

        This method issues requests to the API.

        Args:
            polling_period: How often to poll the API for the task status.

        Returns:
            The final status of the task.
        """
        prev_status = None
        while True:
            status = self.get_status()
            if status != prev_status:
                if status == models.TaskStatusCode.PENDINGINPUT:
                    pass
                elif status == models.TaskStatusCode.SUBMITTED:
                    logging.info(
                        "Task %s successfully queued and waiting to be "
                        "picked-up for execution...", self.id)
                elif status == models.TaskStatusCode.STARTED:
                    logging.info(
                        "Task %s has started and is now running "
                        "remotely.", self.id)
                elif status == models.TaskStatusCode.SUCCESS:
                    logging.info("Task %s completed successfully.", self.id)
                elif status == models.TaskStatusCode.FAILED:
                    logging.info("Task %s failed.". self.id)
                    logging.info("Download the 'stdout.txt' and 'stderr.txt' "
                                 "files with `task.download_outputs()` for "
                                 "more detail.")
                elif status == models.TaskStatusCode.PENDINGKILLED:
                    logging.info("Task %s is being killed.", self.id)
                elif status == models.TaskStatusCode.KILLED:
                    logging.info("Task %s killed.", self.id)
                elif status == models.TaskStatusCode.ZOMBIE:
                    logging.info("The machine was terminated while the task "
                                 "was pending.")
                else:
                    logging.info("An internal error occurred while "
                                 "performing the task.")
            prev_status = status

            if status in _TASK_TERMINAL_STATUSES:
                return status

            time.sleep(polling_period)

    def kill(self) -> None:
        """Kill the task.

        This method issues a request to the API to kill the task. It doesn't
        block waiting for confirmation if the task was killed.
        """
        self._api.kill_task(path_params=self._get_path_params())

    def get_simulator_name(self) -> str:
        # e.g. retrieve openfoam from fvm.openfoam.run_simulation
        return self.get_info()["method_name"].split(".")[1]

    def get_storage_path(self) -> str:
        return self.get_info()["storage_path"]

    def get_output_files_info(self) -> output_contents.OutputContents:
        """Get information of the output files of the task.

        Returns:
            An instance of the OutputContents class, which can be used to
            access info about the output files, such as the archive's size,
            number of files, and information about each file (name, size,
            compressed size). It can also be used to print that information
            in a formatted way.
        """
        api_response = self._api.get_outputs_list(
            path_params=self._get_path_params())

        archive_info = api_response.body

        contents = [
            output_contents.FileInfo(
                name=file_info["name"],
                size=file_info["size"],
                compressed_size=file_info["compressed_size"],
            ) for file_info in archive_info["contents"]
        ]

        return output_contents.OutputContents(
            size=int(archive_info["size"]),
            contents=contents,
        )

    def download_outputs(
        self,
        filenames: Optional[List[str]] = None,
        output_dir: Optional[types.Path] = None,
        uncompress: bool = True,
        rm_downloaded_zip_archive: bool = True,
    ) -> pathlib.Path:
        """Download output files of the task.

        Args:
            filenames: List of filenames to download. If None or empty, all
                files are downloaded.
            output_dir: Directory where to download the files. If None, the
                files are downloaded to the default directory. The default is
                {inductiva.working_dir}/{inductiva.output_dir}/{task_id}.
            uncompress: Whether to uncompress the archive after downloading it.
            rm_downloaded_zip_archive: Whether to remove the archive after
            uncompressing it. If uncompress is False, this argument is ignored.
        """
        api_response = self._api.download_task_output(
            path_params=self._get_path_params(),
            query_params={
                "filename": filenames or [],
            },
            stream=True,
            skip_deserialization=True,
        )
        # use raw urllib3 response instead of the generated client response, to
        # implement our own download logic (with progress bar, first checking
        # the size of the file, etc.)
        response = api_response.response

        if output_dir is None:
            output_dir = files.resolve_path(inductiva.output_dir).joinpath(
                self.id)
        output_dir = files.resolve_path(output_dir)

        if output_dir.exists():
            warnings.warn("Path already exists, files may be overwritten.")
        output_dir.mkdir(parents=True, exist_ok=True)

        zip_path = output_dir.joinpath("output.zip")

        logging.info("Downloading simulation outputs to %s.", zip_path)
        data.download_file(response, zip_path)

        if uncompress:
            logging.info("Uncompressing the outputs to %s.", output_dir)
            data.uncompress_task_outputs(zip_path, output_dir)
            if rm_downloaded_zip_archive:
                zip_path.unlink()
        return output_dir

    class _PathParams(TypedDict):
        """Util class for type checking path params."""
        task_id: str

    def _get_path_params(self) -> _PathParams:
        """Get dictionary with the URL path parameters for API calls."""
        return {"task_id": self.id}

    def get_computation_time(self,
                             fail_if_running: bool = True) -> Optional[float]:
        """Get the time the computation of the task took to complete.

        Returns:
            The time in hh mm ss or None if the task hasn't completed yet.
        """
        info = self.get_info()
        if fail_if_running and self._status not in _TASK_TERMINAL_STATUSES:
            return None
        # start_time may be None if the task was killed before it started
        if info["computation_start_time"] is None:
            return None

        # Format the time to datetime type
        start_time = datetime.datetime.fromisoformat(
            info["computation_start_time"])
        end_time = info.get("computation_end_time")
        if end_time is None:
            end_time = datetime.datetime.now(datetime.timezone.utc)
        else:
            end_time = datetime.datetime.fromisoformat(
                info["computation_end_time"])

        total_seconds = (end_time - start_time).total_seconds()
        return format_utils.seconds_formatter(total_seconds)

    def get_total_time(self, fail_if_running: bool = True) -> Optional[float]:
        """Get the total time the task workflow took to complete.

        Returns:
            The time in hh mm ss or None if the task hasn't completed yet.
        """
        info = self.get_info()
        if fail_if_running and self._status not in _TASK_TERMINAL_STATUSES:
            return None
        # start_time may be None if the task was killed before it started
        if info["input_submit_time"] is None:
            return None

        # Format the time to datetime type
        submitted_time = datetime.datetime.fromisoformat(
            info["input_submit_time"])
        end_time = info.get("end_time")

        if end_time is None:
            end_time = datetime.datetime.now(datetime.timezone.utc)
        else:
            end_time = datetime.datetime.fromisoformat(info["end_time"])

        total_seconds = (end_time - submitted_time).total_seconds()
        return format_utils.seconds_formatter(total_seconds)

    def get_machine_type(self) -> Optional[str]:
        """Get the machine type used in the task.

        Streamlines the process of obtaining the task info, extracting the
        machine type from the comprehensive task info.

        Returns:
            The machine type, or None if a machine hasn't been assigned yet.
        """
        info = self.get_info()
        if info["executer"] is None:
            return None

        machine_info = info["executer"]
        if "vm_type" not in machine_info:
            return "inductiva-machine"

        machine_type = machine_info["vm_type"].split("/")[-1]

        return machine_type

    def get_stdout(self, n_lines: int = 10, verbose: bool = True):
        """Returns tail of stdout.txt file for current task

        Args:
            n_lines: Number of lines to return from the end of the file.
            verbose: Whether to print the contents.
        Returns:
            A list of strings, each string being a line from the stdout."""

        if self.get_status() in (models.TaskStatusCode.PENDINGINPUT,
                                 models.TaskStatusCode.SUBMITTED):
            logging.info("Task has not started yet.")
            return

        api_response = self._api.get_stdout_tail(
            path_params=self._get_path_params(),
            query_params={
                "n_lines": n_lines,
            },
            stream=False,
            skip_deserialization=False,
        )

        if verbose:
            logging.info("Simulation stdout:")
            print("\n")
            for line in api_response.body:
                print(line)

        return api_response.body

    def get_resources_usage(self, n_lines: int = 10, verbose: bool = True):
        """Returns tail of resources_usage.txt file for current task

        Calls the get_resources_tail function. This file is a .csv file
        with each line corresponding to: register_time / memory_usage_percent /
        cpu_usage_percent. The function returns the last n_lines in a list of
        lines.
        Args:
            n_lines: Number of lines to return from the end of the file.
            verbose: Whether to print the contents.
        Returns:
            A list of strings, each string being a line from the stdout."""
        if self.get_status() in (models.TaskStatusCode.PENDINGINPUT,
                                 models.TaskStatusCode.SUBMITTED):
            logging.info("Task did not start yet.")
            return

        api_response = self._api.get_resources_tail(
            path_params=self._get_path_params(),
            query_params={
                "n_lines": n_lines,
            },
            stream=False,
            skip_deserialization=False,
        )

        if verbose:
            logging.info("Current resource usage:")
            print("\n")
            print("Timestamp \t   Memory usage  CPU usage")

            for line in api_response.body:
                date, memory, cpu = line.split(",")

                #Remove the miliseconds from the date
                datetime_date = parser.parse(date)
                truncated_datetime = datetime_date.strftime("%Y-%m-%d %H:%M:%S")

                print(f"{truncated_datetime} \t {float(memory):.3f}\
                       \t {float(cpu):.3f}")

        return api_response.body
