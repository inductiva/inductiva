"""Manage running/completed tasks on the Inductiva API."""
import pathlib
import contextlib
import time
import json
from absl import logging
from typing import Dict, Any, List, Optional, Tuple, Union
from typing_extensions import TypedDict
import datetime
from ..localization import translator as __
import tabulate
import sys

from inductiva import constants
from inductiva.client import exceptions, models
from inductiva import api, types
from inductiva.client.apis.tags import tasks_api
from inductiva.utils import files, format_utils, data
from inductiva.tasks import output_info

import warnings


class Metric:
    """Represents a single metric with a value and a label."""

    def __init__(self, label: str, value=None):
        self.label = label
        self.value = value


class TaskInfo:
    """Represents the task information."""

    class Executer:
        """Encapsulates information about the executer."""

        def __init__(self):
            self.uuid = None
            self.cpu_count_logical = None
            self.cpu_count_physical = None
            self.memory = None
            self.n_mpi_hosts = None
            self.vm_type = None
            self.vm_name = None
            self.host_type = None
            self.error_detail = None

    class TimeMetrics:
        """Encapsulates time-related metrics."""

        def __init__(self):
            self.total_seconds = Metric("Total seconds")
            self.container_image_download_seconds = Metric(
                "Container image download")
            self.queue_time_seconds = Metric("Time in queue")
            self.computation_seconds = Metric("Computation")
            self.input_download_seconds = Metric("Input download")
            self.input_decompression_seconds = Metric("Input decompression")
            self.output_compression_seconds = Metric("Output compression")
            self.output_upload_seconds = Metric("Output upload")

    class DataMetrics:
        """Encapsulates data-related metrics."""

        def __init__(self):
            self.output_zipped_size_bytes = Metric("Size of zipped output")
            self.output_size_bytes = Metric("Size of unzipped output")
            self.output_total_files = Metric("Number of output files")

    def __init__(self, **kwargs):
        """Initialize the instance with attributes from keyword arguments."""
        self._is_running = None
        self._is_terminal = None
        self.task_id = None
        self.status = None
        self.method_name = None
        self.storage_path = None
        self.container_image = None
        self.project = None
        self.create_time = None
        self.input_submit_time = None
        self.start_time = None
        self.computation_start_time = None
        self.computation_end_time = None
        self.end_time = None
        self.time_metrics = self.TimeMetrics()
        self.data_metrics = self.DataMetrics()

        # Set the general attributes
        for key, value in kwargs.items():
            if hasattr(self, key):
                setattr(self, key, value)

        # Set the metrics objects
        metrics_info: dict = kwargs.get("metrics", {})
        self.__update_metrics(self.time_metrics, metrics_info)
        self.__update_metrics(self.data_metrics, metrics_info)

        # Set the executer object
        self.executer = None
        executer_info: dict = kwargs.get("executer", {})
        if executer_info:
            self.executer = self.Executer()
            for key, value in executer_info.items():
                if hasattr(self.executer, key):
                    setattr(self.executer, key, value)

        print("\n", self.__dict__)
        print("\n", self.time_metrics.__dict__)
        print("\n", self.data_metrics.__dict__)
        print("\n", self.executer.__dict__)

    def __update_metrics(
        self,
        metrics_obj: Metric,
        metrics_info: Dict[str, Any],
    ) -> None:
        for key, value in metrics_info.items():
            if key in metrics_obj.__dict__:
                getattr(metrics_obj, key).value = value

    def _format_time_metric(
        self,
        metric_key: str,
        value: Optional[float],
    ) -> str:
        if isinstance(value, float):
            if value >= 60.0:
                value_str = format_utils.seconds_formatter(value)
            else:
                value_str = "{:.2f} s".format(value)

            if metric_key in ("total_seconds",
                              "computation_seconds") and self._is_running:
                value_str += " (Task still running)"

            return value_str

        # Value is None
        if (metric_key == "computation_seconds" and
                self._is_running is False and self._is_terminal is False):
            return "N/A until task is started"

        if metric_key == "container_image_download_seconds":
            return "N/A"

        return value

    @staticmethod
    def _format_data_metric(
        metric_key: str,
        value: Optional[float],
    ) -> Optional[str]:
        if isinstance(value, int) and "bytes" in metric_key:
            return format_utils.bytes_formatter(value)
        return value

    def __str__(self):
        table_format = "plain"
        default_missing_val = "N/A until task ends"

        wall_time_data = [[
            "Wall time:",
            self._format_time_metric(
                "total_seconds",
                self.time_metrics.total_seconds.value,
            ),
        ]]
        wall_time_table = tabulate.tabulate(
            wall_time_data,
            tablefmt=table_format,
        )

        time_metrics_data = [
            [
                f"{metric.label}:",
                self._format_time_metric(metric_key, metric.value),
            ]
            for metric_key, metric in self.time_metrics.__dict__.items()
            if metric_key != "total_seconds"
        ]
        time_metrics_table = tabulate.tabulate(
            time_metrics_data,
            missingval=default_missing_val,
            tablefmt=table_format,
        )

        data_metrics_data = [[
            f"{metric.label}:",
            self._format_data_metric(metric_key, metric.value)
        ] for metric_key, metric in self.data_metrics.__dict__.items()]
        data_metrics_table = tabulate.tabulate(
            data_metrics_data,
            missingval=default_missing_val,
            tablefmt=table_format,
        )

        # Add a tab to the beginning of each line in the tables
        time_metrics_table = "\n".join(
            "\t" + line for line in time_metrics_table.splitlines())
        data_metrics_table = "\n".join(
            "\t" + line for line in data_metrics_table.splitlines())

        table_str = f"\n{wall_time_table}"
        table_str += f"\n\tTime breakdown:\n{time_metrics_table}"
        table_str += f"\nData:\n{data_metrics_table}\n"

        return table_str


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

    FAILED_STATUSES = {
        models.TaskStatusCode.FAILED, models.TaskStatusCode.KILLED,
        models.TaskStatusCode.EXECUTERFAILED,
        models.TaskStatusCode.EXECUTERTERMINATED,
        models.TaskStatusCode.EXECUTERTERMINATEDBYUSER,
        models.TaskStatusCode.SPOTINSTANCEPREEMPTED,
        models.TaskStatusCode.ZOMBIE
    }

    TERMINAL_STATUSES = {models.TaskStatusCode.SUCCESS}.union(FAILED_STATUSES)

    RUNNING_STATUSES = {
        models.TaskStatusCode.PENDINGINPUT, models.TaskStatusCode.STARTED
    }

    KILLABLE_STATUSES = {models.TaskStatusCode.SUBMITTED
                        }.union(RUNNING_STATUSES)

    KILL_VERBOSITY_LEVELS = [0, 1, 2]

    def __init__(self, task_id: str):
        """Initialize the instance from a task ID."""
        self.id = task_id
        self._api = tasks_api.TasksApi(api.get_client())
        self._info = None
        self._status = None

    def is_running(self, cached: bool = False) -> bool:
        """Validate if the task is running.

        This method issues a request to the API.
        """
        return self.get_status(cached) == models.TaskStatusCode.STARTED

    def is_failed(self, cached: bool = False) -> bool:
        """Validate if the task has failed.

        This method issues a request to the API.
        """
        return self.get_status(cached) in self.FAILED_STATUSES

    def is_terminal(self, cached: bool = False) -> bool:
        """Check if the task is in a terminal status.

        This method issues a request to the API.
        """
        return self.get_status(cached) in self.TERMINAL_STATUSES

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

    def get_status(self, cached: bool = False) -> models.TaskStatusCode:
        """Get status of the task.

        This method issues a request to the API.
        """
        # If the task is in a terminal status and we already have the status,
        # return it without refreshing it from the API.
        # Can't call is_terminal because it calls get status (infinite loop)
        if (self._status is not None and
                self._status in self.TERMINAL_STATUSES) or cached is True:
            return self._status

        resp = self._api.get_task_status(self._get_path_params())

        return models.TaskStatusCode(resp.body["status"])

    @property
    def info(self) -> TaskInfo:
        """Get information about the task.

        It contains cached information about the task from the latest call to
        `get_info`, therefore it can be outdated.
        """
        return self._info

    def get_info(self) -> TaskInfo:
        """Get a dictionary with information about the task.

        Information includes e.g., "task_id", "status", timestamps
        ("create_time", "input_submit_time, "start_time", "end_time"),
        among others.

        This method issues a request to the API.
        """
        params = self._get_path_params()
        resp = self._api.get_task(params, skip_deserialization=True).response

        info = json.loads(resp.data.decode("utf-8"))
        status = models.TaskStatusCode(info["status"])

        self._info = TaskInfo(**info)
        self._status = status

        return self._info

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
                    logging.info("Task %s failed.", self.id)
                    logging.info("Download the 'stdout.txt' and 'stderr.txt' "
                                 "files with `task.download_outputs()` for "
                                 "more detail.")
                elif status == models.TaskStatusCode.PENDINGKILL:
                    logging.info("Task %s is being killed.", self.id)
                elif status == models.TaskStatusCode.KILLED:
                    logging.info("Task %s killed.", self.id)
                elif status == models.TaskStatusCode.ZOMBIE:
                    logging.info("The machine was terminated while the task "
                                 "was pending.")
                elif status == models.TaskStatusCode.EXECUTERFAILED:
                    info = self.get_info()
                    detail = info.executer.error_detail
                    logging.info("The remote process running the task failed:")
                    if detail:
                        logging.info(" > Message: %s", detail)
                    else:
                        logging.info(" > No error message available.")

                else:
                    logging.info(
                        "An internal error occurred with status %s "
                        "while performing the task.", status)
            prev_status = status

            if self.is_terminal():
                return status

            time.sleep(polling_period)

    def _send_kill_request(self, max_api_requests: int) -> None:
        """Send a kill request to the API.
        If the api request fails, it will retry until
        max_api_requests is 0 raising a RuntimeError.
        Args:
            max_api_requests (int): maximum number of api requests to send
        """
        while max_api_requests > 0:
            max_api_requests -= 1
            try:
                if self.is_terminal():
                    break

                path_params = self._get_path_params()
                self._api.kill_task(path_params=path_params)
                break
            except exceptions.ApiException as exc:
                if max_api_requests == 0:
                    raise RuntimeError(
                        "Something went wrong while sending"
                        " the kill command. Please try again later.") from exc
                time.sleep(constants.TASK_KILL_RETRY_SLEEP_SEC)

    def _check_if_pending_kill(
            self,
            wait_timeout: Union[float,
                                int]) -> Tuple[bool, models.TaskStatusCode]:
        """Check if the task is in the PENDINGKILL state.
        This method keeps checking the status of the task until it is no longer
        in the PENDINGKILL state or until the timeout is reached.
        Args:
            wait_timeout (int, float): number of seconds to wait for the
            state to leave PENDINGKILL.
        Returns:
            A tuple with a boolean indicating whether the timeout was reached
            and the status of the task.
        """
        success = True
        start = time.time()

        while (status :=
               self.get_status()) == models.TaskStatusCode.PENDINGKILL:
            if (time.time() - start) > wait_timeout:
                success = False
                break

            time.sleep(constants.TASK_KILL_RETRY_SLEEP_SEC)
        return success, status

    def kill(self,
             wait_timeout: Optional[Union[float, int]] = None,
             verbosity_level: int = 2) -> Union[bool, None]:
        """Request a task to be killed.

        This method requests that the current task is remotely killed.
        If `wait_timeout` is None (default), the kill request is sent to the
        backend and the method returns. However, if `wait_timeout` is
        a positive number, the method waits up to `wait_timeout` seconds
        to ensure that the task transitions to the KILLED state.
        Args:
            wait_timeout (int, float): Optional - number of seconds to wait
            for the kill command or None if only the request is to be sent.
            verbosity_level (int): Optional. the verbosity of the logs when the
            task signal is sent and when the task is killed. Verbosity 0
            produces no outputs, 1 produces minimal outputs, and 2 (Default)
            produces extensive outputs.
        Returns:
            - None if `wait_timeout` is None and the kill request was
              successfully sent;
            - True if `wait_timeout`> 0 and the task successfully transitioned
              to the `KILLED` state within `wait_timeout` seconds;
            - False if `wait_timeout` > 0 but the task didn't transition
              to the `KILLED` state within `wait_timeout` seconds;

        """
        if wait_timeout is not None:
            if not isinstance(wait_timeout, (float, int)):
                raise TypeError("Wait timeout must be a number.")
            if wait_timeout <= 0.0:
                raise ValueError("Wait timeout must be a positive number.")

        if verbosity_level not in self.KILL_VERBOSITY_LEVELS:
            raise ValueError(f"Verbosity {verbosity_level} level not allowed. "
                             f"Choose from {self.KILL_VERBOSITY_LEVELS}")

        self._send_kill_request(constants.TASK_KILL_MAX_API_REQUESTS)

        if wait_timeout is None:
            logging.info(
                __("task-kill-request-sent" + f"-{verbosity_level}", self.id))
            return None

        success, status = self._check_if_pending_kill(wait_timeout)

        if status != models.TaskStatusCode.KILLED:
            success = False
            if status == models.TaskStatusCode.PENDINGKILL:
                logging.error(
                    "Unable to ensure that task %s transitioned to the KILLED "
                    "state after %f seconds. The status of the task is %s.",
                    self.id, wait_timeout, status)
            else:
                logging.error(
                    "Task is already in a terminal state and cannot be killed. "
                    "Current task status is %s.", status)

        if success:
            if verbosity_level == 2:
                logging.info("Successfully killed task %s.", self.id)
            elif verbosity_level == 1:
                logging.info("%s killed.")

        return success

    def get_simulator_name(self) -> str:
        # e.g. retrieve openfoam from fvm.openfoam.run_simulation
        return self.get_info().method_name.split(".")[1]

    def get_storage_path(self) -> str:
        return self.get_info().storage_path

    def get_output_info(self) -> output_info.TaskOutputInfo:
        """Get information about the output files of the task.

        Returns:
            An instance of the OutputInfo class, which can be used to
            access info about the output archive (number of files, total
            compressed size, total uncompressed size) and information about
            each file (name, size, compressed size). It can also be used to
            print that information in a formatted way.
        """
        api_response = self._api.get_outputs_list(
            path_params=self._get_path_params())

        archive_info = api_response.body

        output_files = [
            output_info.FileInfo(
                name=file_info["name"],
                size=int(file_info["size"]),
                compressed_size=int(file_info["compressed_size"]),
            ) for file_info in archive_info["contents"]
        ]

        return output_info.TaskOutputInfo(
            total_size_bytes=int(archive_info["size"]),
            files=output_files,
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
                {inductiva.get_output_dir()}/{output_dir}/{task_id}.
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
            output_dir = files.resolve_output_path(self.id)
        else:
            output_dir = files.resolve_output_path(output_dir)

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

    def _get_duration(
        self,
        start_attribute: str,
        metric_attribute: str,
        cached: bool,
    ) -> Optional[float]:
        """Get the duration of a task phase.

        Args:
            start_attribute: The attribute containing the start time.
            metric_attribute: The attribute containing the duration if the task
                has ended.
            cached: Whether to use the cached info or fetch the latest info.

        Returns:
            The duration in seconds, or None if the start or end time is None.
        """
        info: TaskInfo = self.get_info() if not cached else self.info

        # The task has ended and the metric is available
        metric = getattr(info.time_metrics, metric_attribute)
        if metric.value is not None:
            return metric.value

        # The task has ended but the metric is not available
        if self.is_terminal():
            return None

        # The task is still running
        start_time = getattr(info, start_attribute)

        # start time may be None if the task was killed before it started
        if start_time is None:
            return None

        # Format the time to datetime type
        start_time = datetime.datetime.fromisoformat(start_time)
        end_time = datetime.datetime.now(datetime.timezone.utc)

        return (end_time - start_time).total_seconds()

    def get_computation_time(self, cached: bool = False) -> Optional[float]:
        """Get the time the computation of the task took to complete.

        Returns:
            The task computation time if the task is already started or in a
            terminal state, or None otherwise.
        """
        return self._get_duration(
            start_attribute="computation_start_time",
            metric_attribute="computation_seconds",
            cached=cached,
        )

    def get_total_time(self, cached: bool = False) -> Optional[float]:
        """Get the total time the task workflow took to complete.

        Returns:
            The task total duration since it was created, or None if the
            metric is not available or can't be computed.
        """
        return self._get_duration(
            start_attribute="create_time",
            metric_attribute="total_seconds",
            cached=cached,
        )

    def get_machine_type(self) -> Optional[str]:
        """Get the machine type used in the task.

        Streamlines the process of obtaining the task info, extracting the
        machine type from the comprehensive task info.

        Returns:
            The machine type, or None if a machine hasn't been assigned yet.
        """
        info: TaskInfo = self.get_info()
        if info.executer is None:
            return None

        machine_provider = info.executer.host_type
        machine_type = info.executer.vm_type.split("/")[-1]

        return machine_provider + "-" + machine_type

    def print_summary(self, fhandle=sys.stdout):
        info: TaskInfo = self.get_info()

        # Update running info
        info._is_running = self.is_running(cached=True)
        info._is_terminal = self.is_terminal(cached=True)

        # Update the duration metrics if the task is still running
        if info._is_running is True:
            info.time_metrics.total_seconds.value = self.get_total_time(
                cached=True)
            info.time_metrics.computation_seconds.value = \
                self.get_computation_time(cached=True)

        print(info, file=fhandle)
