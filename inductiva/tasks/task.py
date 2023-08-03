"""Manage running/completed tasks on the Inductiva API."""
import pathlib
import shutil
import time
from absl import logging
from typing import Dict, Any, List, Optional
from typing_extensions import TypedDict
import datetime
import inductiva
from inductiva.client.models import TaskStatusCode
from inductiva import api
from inductiva.client.apis.tags.tasks_api import TasksApi
from inductiva.utils import files
from inductiva.utils import data
from inductiva.utils import output_contents


class Task:
    """Represents a running/completed task on the Inductiva API.

    Example usage:
        task = scenario.simulate_async(...)
        final_status = task.wait()
        info = task.get_info() # dictionary with info about the task
        task.download_output() # download the output of the task

    Attributes:
        id: The task ID.
        _api: Instance of TasksApi (from the generated API client).
    """

    def __init__(self, task_id: str):
        """Initialize the instance from a task ID."""
        self.id = task_id
        self._api = TasksApi(api.get_client())

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

    def __exit__(self, exc_type, exc_value, traceback):
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

    def get_status(self) -> TaskStatusCode:
        """Get status of the task.

        This method issues a request to the API.
        """
        resp = self._api.get_task_status(self._get_path_params())
        return TaskStatusCode(resp.body["status"])

    def get_info(self) -> Dict[str, Any]:
        """Get a dictionary with information about the task.

        Information includes e.g., "task_id", "status", timestamps
        ("create_time", "input_submit_time, "start_time", "end_time"),
        among others.

        This method issues a request to the API.
        """
        params = self._get_path_params()
        resp = self._api.get_task(params, skip_deserialization=True).response
        return resp.json()

    def wait(self, polling_period: int = 5) -> TaskStatusCode:
        """Wait for the task to complete.

        This method issues requests to the API.

        Args:
            polling_period: How often to poll the API for the task status.

        Returns:
            The final status of the task.
        """
        terminal_statuses = {
            TaskStatusCode.SUCCESS,
            TaskStatusCode.FAILED,
            TaskStatusCode.KILLED,
            TaskStatusCode.EXECUTERFAILED,
            TaskStatusCode.EXECUTERTERMINATED,
            TaskStatusCode.SPOTINSTANCEPREEMPTED,
        }

        prev_status = None
        while True:
            status = self.get_status()
            if status != prev_status:
                if status == TaskStatusCode.PENDINGINPUT:
                    pass
                elif status == TaskStatusCode.SUBMITTED:
                    logging.info("Waiting for resources...")
                elif status == TaskStatusCode.STARTED:
                    logging.info("The task is being executed remotely.")
                elif status == TaskStatusCode.SUCCESS:
                    logging.info("Task completed successfully.")
                elif status == TaskStatusCode.FAILED:
                    logging.info("Task failed.")
                    logging.info("Download the task output and check the "
                                 "'stdout.txt' and 'stderr.txt' files for "
                                 "more information.")
                elif status == TaskStatusCode.KILLED:
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

    def get_outputs_info(self):
        api_response = self._api.get_outputs_list(
            path_params=self._get_path_params(),)

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
        output_dir: Optional[pathlib.Path] = None,
        uncompress: bool = True,
        rm_archive: bool = True,
    ) -> pathlib.Path:
        api_response = self._api.download_task_output(
            path_params=self._get_path_params(),
            query_params={
                "filename": filenames or [],
            },
            stream=True,
            skip_deserialization=True,
        )
        # use raw urllib3 response instead of generated client response, to
        # implement our own download logic (with progress bar, first checking
        # the size of the file, etc.)
        response = api_response.response

        if output_dir is None:
            output_dir = files.resolve_path(inductiva.output_dir).joinpath(
                self.id)

        if output_dir.exists():
            shutil.rmtree(output_dir)

        output_dir.mkdir(parents=True)

        zip_path = output_dir.joinpath("output.zip")

        data.download_file(response, zip_path)

        response.release_conn()

        if uncompress:
            data.uncompress_task_outputs(zip_path, output_dir)
            if rm_archive:
                zip_path.unlink()

        return output_dir

    class _PathParams(TypedDict):
        """Util class for type checking path params."""
        task_id: str

    def _get_path_params(self) -> _PathParams:
        "Get dictionary representing the URL path parameters for API calls."
        return {"task_id": self.id}

    def get_execution_time(self) -> float:
        """Get the time the task took to complete.

        Returns:
            The time in seconds.
        """

        if self.get_status() != TaskStatusCode.SUCCESS:
            raise RuntimeError("Task is not completed.")

        params = self._get_path_params()
        info = dict(self._api.get_task(params).body)

        #Format the time to datetime type
        end_time = datetime.datetime.strptime(str(info["end_time"]),
                                              "%Y-%m-%dT%H:%M:%S.%f+00:00")
        start_time = datetime.datetime.strptime(str(info["start_time"]),
                                                "%Y-%m-%dT%H:%M:%S.%f+00:00")

        return (end_time - start_time).total_seconds()
