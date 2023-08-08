"""Manage running/completed tasks on the Inductiva API."""
import pathlib
import time
import json
from absl import logging
from typing import Dict, Any, List, Optional
from typing_extensions import TypedDict
import datetime
import inductiva
from inductiva.client import models
from inductiva import api
from inductiva.client.apis.tags import tasks_api
from inductiva.utils import files
from inductiva.utils import data
from inductiva.utils import output_contents
from inductiva import types
import warnings


class Task:
    """Represents a running/completed task on the Inductiva API.

    Example usage:
        task = scenario.simulate_async(...)
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
        self._output_class = None

    def __enter__(self):
        """Enter context manager for managing a blocking execution.

        If an exception/ctrl+c is caught while in the context manager, the task
        is killed.

        Usage:
            task = scenario.simulate_async(...)
            with task:
                # If an exception happens here or ctrl+c is pressed, the task
                # will be killed.
                task.wait()

        """
        return self

    def __exit__(self, exc_type, exc_value, traceback) -> bool:
        """Exit context manager killing the task if an exception was raised."""
        del traceback  # unused
        del exc_type  # unused

        if exc_value is None:
            return True

        if isinstance(exc_value, KeyboardInterrupt):
            logging.info("Caught SIGINT: terminating blocking task...")
        elif exc_value is not None:
            logging.info("Caught exception: terminating blocking task...")

        self.kill()
        return False

    def get_status(self) -> models.TaskStatusCode:
        """Get status of the task.

        This method issues a request to the API.
        """
        resp = self._api.get_task_status(self._get_path_params())
        return models.TaskStatusCode(resp.body["status"])

    def get_info(self) -> Dict[str, Any]:
        """Get a dictionary with information about the task.

        Information includes e.g., "task_id", "status", timestamps
        ("create_time", "input_submit_time, "start_time", "end_time"),
        among others.

        This method issues a request to the API.
        """
        params = self._get_path_params()
        resp = self._api.get_task(params, skip_deserialization=True).response

        return json.loads(resp.data.decode("utf-8"))

    def wait(self, polling_period: int = 5) -> models.TaskStatusCode:
        """Wait for the task to complete.

        This method issues requests to the API.

        Args:
            polling_period: How often to poll the API for the task status.

        Returns:
            The final status of the task.
        """
        terminal_statuses = {
            models.TaskStatusCode.SUCCESS,
            models.TaskStatusCode.FAILED,
            models.TaskStatusCode.KILLED,
            models.TaskStatusCode.EXECUTERFAILED,
            models.TaskStatusCode.EXECUTERTERMINATED,
            models.TaskStatusCode.SPOTINSTANCEPREEMPTED,
        }

        prev_status = None
        while True:
            status = self.get_status()
            if status != prev_status:
                if status == models.TaskStatusCode.PENDINGINPUT:
                    pass
                elif status == models.TaskStatusCode.SUBMITTED:
                    logging.info("Waiting for resources...")
                elif status == models.TaskStatusCode.STARTED:
                    logging.info("The task is being executed remotely.")
                elif status == models.TaskStatusCode.SUCCESS:
                    logging.info("Task completed successfully.")
                elif status == models.TaskStatusCode.FAILED:
                    logging.info("Task failed.")
                    logging.info("Download the task output and check the "
                                 "'stdout.txt' and 'stderr.txt' files for "
                                 "more information.")
                elif status == models.TaskStatusCode.KILLED:
                    logging.info("Task killed.")
                else:
                    logging.info(
                        "An internal error occurred while performing the task.")
            prev_status = status

            if status in terminal_statuses:
                return status

            time.sleep(polling_period)

    def kill(self) -> None:
        """Kill the task.

        This method issues a request to the API to kill the task. It doesn't
        block waiting for confirmation if the task was killed.
        """
        self._api.kill_task(path_params=self._get_path_params())

    def set_output_class(self, output_class):
        """Set the output class of the task."""
        self._output_class = output_class

    def get_output(self,
                   output_dir: Optional[types.Path] = None) -> pathlib.Path:
        """Get the output files produced by the task. If the task 
        is not finished the function will wait until it is to download
        the output folder contents. 
        Returns:
            The output path of the task.
        Example:
            task = Task("task_id") 
            output_path = task.get_output()  
            # prints:
            100%|██████████| 1.64G/1.64G [00:32<00:00, 55.1MiB/s]   
        """
        self.wait()
        output_dir = self.download_outputs(output_dir=output_dir)

        if self._output_class:
            return self._output_class(output_dir)

        return output_dir

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
        rm_archive: bool = True,
    ) -> pathlib.Path:
        """Download output files of the task.

        Args:
            filenames: List of filenames to download. If None or empty, all
                files are downloaded.
            output_dir: Directory where to download the files. If None, the
                files are downloaded to the default directory. The default is
                {inductiva.working_dir}/{inductiva.output_dir}/{task_id}.
            uncompress: Whether to uncompress the archive after downloading it.
            rm_archive: Whether to remove the archive after uncompressing it.
                If uncompress is False, this argument is ignored.
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

        data.download_file(response, zip_path)

        if uncompress:
            data.uncompress_task_outputs(zip_path, output_dir)
            if rm_archive:
                zip_path.unlink()

        return output_dir

    class _PathParams(TypedDict):
        """Util class for type checking path params."""
        task_id: str

    def _get_path_params(self) -> _PathParams:
        """Get dictionary with the URL path parameters for API calls."""
        return {"task_id": self.id}

    def get_execution_time(self) -> float:
        """Get the time the task took to complete.

        Returns:
            The time in seconds.
        """

        if self.get_status() != models.TaskStatusCode.SUCCESS:
            raise RuntimeError("Task is not completed.")

        params = self._get_path_params()
        info = dict(self._api.get_task(params).body)

        #Format the time to datetime type
        end_time = datetime.datetime.strptime(str(info["end_time"]),
                                              "%Y-%m-%dT%H:%M:%S.%f+00:00")
        start_time = datetime.datetime.strptime(str(info["start_time"]),
                                                "%Y-%m-%dT%H:%M:%S.%f+00:00")

        return (end_time - start_time).total_seconds()

    def get_machine_type(self) -> str:
        """Get the machine type used in the task.

        Streamlines the process of obtaining the task info, extracting the
        machine type from the comprehensive task info.

        Returns:
            The machine type.
        """

        params = self._get_path_params()
        task_info = dict(self._api.get_task(params).body)

        machine_info = dict(task_info["executer"])
        machine_type = machine_info["vm_type"].split("/")[-1]

        return machine_type
