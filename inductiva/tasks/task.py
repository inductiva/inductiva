"""Manage running/completed tasks on the Inductiva API."""
import io
import sys
import time
import json
import asyncio
import pathlib
import logging
import datetime
import contextlib
import nest_asyncio
from typing_extensions import TypedDict
from typing import AsyncGenerator, Callable, Dict, Any, List, Optional, TextIO, Tuple, Union

import urllib3
import tabulate
from dataclasses import dataclass

from ..localization import translator as __

import inductiva
from inductiva import storage
from inductiva import constants
from inductiva.client import exceptions, models
from inductiva import api
from inductiva.client.apis.tags import tasks_api
from inductiva.utils import files, format_utils, data
from inductiva.tasks import output_info
from inductiva.tasks.file_tracker import Operations, FileTracker

import warnings


@dataclass
class Metric:
    """Represents a single metric with a value and a label.
    
    :meta private:
    """
    label: str
    value: Optional[float] = None


class TaskInfo:
    """Represents the task information.

    :meta private:
    """

    MISSING_UNTIL_TASK_STARTED = "N/A until task is started"
    MISSING_UNTIL_TASK_ENDED = "N/A until task ends"

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
            self.input_upload_seconds = Metric("Input upload")
            self.queue_time_seconds = Metric("Time in queue")
            self.container_image_download_seconds = Metric(
                "Container image download")
            self.input_download_seconds = Metric("Input download")
            self.input_decompression_seconds = Metric("Input decompression")
            self.computation_seconds = Metric("Computation")
            self.output_upload_seconds = Metric("Output upload")

    class DataMetrics:
        """Encapsulates data-related metrics."""

        def __init__(self):
            self.output_zipped_size_bytes = Metric("Size of zipped output")
            self.output_size_bytes = Metric("Size of unzipped output")
            self.output_total_files = Metric("Number of output files")

    def __init__(self, **kwargs):
        """Initialize the instance with attributes from keyword arguments."""
        self.is_submitted = None
        self.is_running = None
        self.is_terminal = None
        self.task_id = None
        self.status = None
        self.status_alias = None
        self.simulator = None
        self.storage_path = None
        self.storage_input_path = None
        self.storage_output_path = None
        self.container_image = None
        self.project = None
        self.create_time = None
        self.input_submit_time = None
        self.start_time = None
        self.computation_start_time = None
        self.computation_end_time = None
        self.end_time = None
        self.estimated_computation_cost = None
        self.time_metrics = self.TimeMetrics()
        self.data_metrics = self.DataMetrics()
        self.status_history = []
        self._kwargs = kwargs

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

        # Update running info
        self.is_submitted = self.status == models.TaskStatusCode.SUBMITTED
        self.is_running = self.status in (
            models.TaskStatusCode.STARTED,
            models.TaskStatusCode.COMPUTATIONSTARTED,
            models.TaskStatusCode.COMPUTATIONENDED)
        self.is_terminal = kwargs.get("is_terminated", False)

    def to_dict(self) -> Dict[str, Any]:
        """Return a dictionary with the task information."""
        return self._kwargs

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
            value_str = f"{value:.2f} s"

            if metric_key == "computation_seconds" and self.is_running:
                value_str += " (Task still running)"

            return value_str

        # Value is None if it is not float
        if self.is_terminal:
            if metric_key == "container_image_download_seconds":
                # If the container image is already present in the local cache
                # the download is skipped, therefore the metric does not exist
                return "N/A (used cached image)"
            if metric_key == "computation_seconds":
                # The task might have ended but the metric is not available in
                # the database yet
                return "N/A"

        if metric_key == "container_image_download_seconds":
            # If the task has not ended the local cache image may be used or
            # the download could be in progress
            return "N/A"

        if metric_key in ("output_upload_seconds",):
            return self.MISSING_UNTIL_TASK_ENDED

        # None values will be replaced with a default missing value message
        return value

    @staticmethod
    def _format_data_metric(
        metric_key: str,
        value: Optional[float],
    ) -> Optional[str]:
        if isinstance(value, int) and "bytes" in metric_key:
            return format_utils.bytes_formatter(value)
        return value

    def __repr__(self) -> str:
        return str(self)

    def get_task_time(self) -> str:
        """Computes the task time.

        Computes the task time based on the creation and end time.
        """
        #has not started yet
        if self.create_time is None:
            task_time = "Task not started yet (N/A)"
        #still running
        elif self.end_time is None:
            now = datetime.datetime.now(datetime.timezone.utc)
            create_time = datetime.datetime.fromisoformat(self.create_time)
            task_time = now - create_time
            task_time = self._format_time_metric("total_seconds",
                                                 task_time.total_seconds())
        else:
            create_time = datetime.datetime.fromisoformat(self.create_time)
            end_time = datetime.datetime.fromisoformat(self.end_time)
            task_time = end_time - create_time
            task_time = self._format_time_metric("total_seconds",
                                                 task_time.total_seconds())
        #if task has ended we use the metric
        return task_time

    def __str__(self):
        table_format = "plain"

        data_metrics_data = [[
            f"{metric.label}:",
            self._format_data_metric(metric_key, metric.value)
        ] for metric_key, metric in self.data_metrics.__dict__.items()]
        data_metrics_table = tabulate.tabulate(
            data_metrics_data,
            missingval=self.MISSING_UNTIL_TASK_ENDED,
            tablefmt=table_format,
        )

        data_metrics_table = "\n".join(
            "\t" + line for line in data_metrics_table.splitlines())

        table_str = f"\nTask status: {self.status_alias}\n"
        if self.executer and self.executer.error_detail:
            table_str += f"\n\tStatus detail: {self.executer.error_detail}"

        table_str += "\nTimeline:\n"
        for item in self.status_history:
            formatted_timestamp = format_utils.datetime_formatter(
                item["timestamp"])

            if item["end_timestamp"]:
                start_time = datetime.datetime.fromisoformat(item["timestamp"])
                end_time = datetime.datetime.fromisoformat(
                    item["end_timestamp"])
                duration_time = end_time - start_time
                duration = f"{round(duration_time.total_seconds(), 3)} s"
            elif not self.is_terminal:
                duration = "(ongoing)"
            else:
                duration = ""

            status = item["alias"]
            table_str += (f"\t{status:<25} at "
                          f"{formatted_timestamp:<20} {duration}\n")

            for index, sub_item in enumerate(item.get("operations", [])):

                if index + 1 == len(item.get("operations", [])):
                    ascii_char = "└"
                else:
                    ascii_char = "├"

                start_time = datetime.datetime.fromisoformat(
                    sub_item["start_timestamp"])
                end_time = datetime.datetime.fromisoformat(
                    sub_item["end_timestamp"]
                ) if sub_item["end_timestamp"] else datetime.datetime.now(
                    datetime.timezone.utc)
                duration_time = end_time - start_time
                duration = f"{round(duration_time.total_seconds(), 3)} s"

                if ("attributes" in sub_item and
                        "command" in sub_item["attributes"]):
                    table_str += (f"\t\t{ascii_char}> {duration:<15} "
                                  f"{sub_item['attributes']['command']}\n")

        table_str += f"\nData:\n{data_metrics_table}\n"
        if self.estimated_computation_cost:
            estimated_cost = format_utils.currency_formatter(
                self.estimated_computation_cost,)
        else:
            estimated_cost = "N/A"
        table_str += ("\nEstimated computation cost (US$): "
                      f"{estimated_cost}\n\n")
        table_str += ("Go to "
                      f"https://console.inductiva.ai/tasks/{self.task_id} "
                      "for more details.")

        return table_str


class Task:
    """Represents a running/completed task on the Inductiva API.

    Example usage:
    
    .. code-block:: python

        task = scenario.simulate(...)
        final_status = task.wait()
        info = task.get_info() # dictionary with info about the task
        task.download_outputs(
            filenames=["file1.txt", "file2.dat"] # download only these files
        )
    """

    FAILED_STATUSES = {
        models.TaskStatusCode.FAILED,
        models.TaskStatusCode.KILLED,
        models.TaskStatusCode.EXECUTERFAILED,
        models.TaskStatusCode.EXECUTERTERMINATED,
        models.TaskStatusCode.EXECUTERTERMINATEDBYUSER,
        models.TaskStatusCode.SPOTINSTANCEPREEMPTED,
        models.TaskStatusCode.ZOMBIE,
        models.TaskStatusCode.EXECUTERTERMINATEDTTLEXCEEDED,
        models.TaskStatusCode.TTLEXCEEDED,
    }

    RUNNING_STATUSES = {
        models.TaskStatusCode.PENDINGINPUT, models.TaskStatusCode.STARTED,
        models.TaskStatusCode.COMPUTATIONSTARTED,
        models.TaskStatusCode.COMPUTATIONENDED
    }

    KILLABLE_STATUSES = {models.TaskStatusCode.SUBMITTED
                        }.union(RUNNING_STATUSES)

    KILL_VERBOSITY_LEVELS = [0, 1, 2]

    STANDARD_OUTPUT_FILES = ["stdout.txt", "stderr.txt"]

    def __init__(self, task_id: str):
        """Initialize the instance from a task ID."""
        self.id = task_id
        self._api = tasks_api.TasksApi(api.get_client())
        self.file_tracker = FileTracker()
        self._info = None
        self._status = None
        self._tasks_ahead: Optional[int] = None
        self._summary = None
        # Internal state to track if the method was called from the wait method
        self._called_from_wait = False

    def is_running(self) -> bool:
        """Validate if the task is running.

        This method issues a request to the API.
        """
        return self.get_status() in (models.TaskStatusCode.STARTED,
                                     models.TaskStatusCode.COMPUTATIONSTARTED,
                                     models.TaskStatusCode.COMPUTATIONENDED)

    def is_failed(self) -> bool:
        """Validate if the task has failed.

        This method issues a request to the API.
        """
        return self.get_status() in self.FAILED_STATUSES

    def is_terminal(self) -> bool:
        """Check if the task is in a terminal status.

        This method issues a request to the API.
        """
        self.get_info()
        return self.info.is_terminal

    @classmethod
    def from_api_info(cls, info: Dict[str, Any]) -> "Task":

        task = cls(info["task_id"])
        task._info = TaskInfo(**info)
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

        This method issues a request to the API and updates the task info
        is_terminal. The api call to get the status now returns the status and
        is_terminated.
        """
        # If the task is in a terminal status and we already have the status,
        # return it without refreshing it from the API.
        if (self._status is not None and self._info.is_terminal):
            return self._status

        resp = self._api.get_task_status(self._get_path_params())

        status = models.TaskStatusCode(resp.body["status"])
        self._status = status

        #updates the info.is_terminal when getting the status
        self._info.is_terminal = resp.body.get(
            "is_terminated",
            self.info.is_terminal,
        )

        queue_position = resp.body.get("position_in_queue", None)
        if queue_position is not None:
            self._tasks_ahead = queue_position.get("tasks_ahead", None)

        return status

    def get_position_in_queue(self) -> Optional[int]:
        """Get the position of the task in the queue.

        This method issues a request to the API.
        """
        _ = self.get_status()

        return self._tasks_ahead

    def _get_last_n_lines_from_file(self, file_path: pathlib.Path,
                                    n: int) -> List[str]:
        """Gets the last n lines from a file.

        This method returns a list with the last n lines from a file.
        Args:
            file_path: The path to the file.
            n: The number of lines to return.
        Returns:
            A list with the last n lines from the file.
        """
        with open(file_path, "r", encoding="utf-8") as file:
            lines = file.readlines()
        return lines[-n:]

    @property
    def info(self) -> TaskInfo:
        """Get information about the task.

        It contains cached information about the task from the latest call to
        `get_info`, therefore it can be outdated.
        """
        if self._info is None:
            return self.get_info()
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

    def _setup_queue_message(self, is_tty: bool, duration: str) -> str:
        if self._tasks_ahead == 0:
            s = f"Task {self.id} is about to start. {duration}"
        else:
            s = (f"Tasks ahead in queue: "
                 f"{self._tasks_ahead} {duration}")
        if not is_tty:
            # We do this because notebooks do not support some escape sequences
            # like the one used to clear the line. So we need to move the cursor
            # to the beginning of the line and overwrite the previous message.
            max_line_length = 73
            return s.ljust(max_line_length, " ")
        return s

    def _format_list_of_lines(self,
                              lines: List[str],
                              file: str,
                              endl: Optional[str] = "",
                              header: Optional[bool] = True) -> str:
        """Formats a list of lines with color.

        This method formats a list of lines with a color and adds a header and
        footer to the list. The color is used to differentiate between stdout
        and stderr.

        Example:
        ┌ (last 10 lines from stderr)
        │ #9  0x7f4ce0941c2d in ???
        │ #10  0x7f4ceb033199 in ???
        │ #11  0x7f4cec45a864 in ???
        │ #12  0x7f4cec4913a6 in ???
        │ #13  0x7f4cecceb940 in ???
        │ #14  0x4084ae in ???
        └
        Args:
            lines: A list of strings to format.
            file: The name of the file. Must be "stdout.txt" or "stderr.txt".
        """

        color_code = "\033[31m" if file == "stderr.txt" else "\033[34m"
        reset_color = "\033[0m"

        if not inductiva.ansi_enabled:
            color_code = ""
            reset_color = ""

        n = len(lines)

        new_lst = [f"{color_code}│{reset_color}{line}{endl}" for line in lines]

        if header:
            new_lst.insert(
                0, f"{color_code}┌ (last {n} lines from {file}){reset_color}\n")
            new_lst.append(f"{color_code}└{reset_color}\n")

        return "".join(new_lst)

    def _format_directory_listing(self,
                                  directories,
                                  current_path="",
                                  indent=0) -> str:
        """Formats a dictionary with directory information.

        This method formats a dictionary with directory information and
        returns a string with the formatted data.

        Args:
            directories: A dictionary with directory information.
            current_path: The current path of the search.
            indent: The current indentation level.
        """

        def build_prefix(indent, is_last):
            """Creates the correct prefix for the current level."""
            color_code = "\033[34m" if inductiva.ansi_enabled else ""
            reset_color = "\033[0m" if inductiva.ansi_enabled else ""
            if indent == 0:
                return ""
            branch = "└── " if is_last else "├── "
            return color_code + "│   " * (indent - 1) + branch + reset_color

        contents = ""
        for index, item in enumerate(directories):
            is_last_item = index == len(directories) - 1
            prefix = build_prefix(indent, is_last_item)

            if isinstance(item, dict):
                for dir_name, dir_contents in item.items():
                    folder_path = f"{current_path}/{dir_name}".strip("/")
                    contents += f"{prefix}{folder_path}/\n"
                    contents += self._format_directory_listing(
                        dir_contents, folder_path, indent + 1)
            else:
                file_path = f"{current_path}/{item}".strip("/")
                contents += f"{prefix}{file_path}\n"

        return contents

    def _print_failed_message(self, out_dir: str) -> None:
        """Prints the messages when a task fails.

        This method prints the last N lines of the stdout and stderr files
        and logs a message to the user to inspect the files for more info.
        Args:
            out_dir: The directory where the files are stored.
        """
        n = constants.TASK_FAILED_LINES_TO_DUMP
        std_out_lines = self._get_last_n_lines_from_file(
            f"{out_dir}/stdout.txt", n)
        std_err_lines = self._get_last_n_lines_from_file(
            f"{out_dir}/stderr.txt", n)

        logging.error("")

        formatted = self._format_list_of_lines(std_out_lines, "stdout.txt")
        logging.error(formatted)
        formatted = self._format_list_of_lines(std_err_lines, "stderr.txt")
        logging.error(formatted)

        logging.error(
            "Please inspect the stdout.txt and stderr.txt files at %s\n"
            "For more information.", out_dir)

    def _handle_status_change(self, status: models.TaskStatusCode,
                              description: str) -> None:
        """Handle a status change.

        Prints a message to the user when the status of the task changes.
        """

        if status == models.TaskStatusCode.EXECUTERFAILED:
            info = self.get_info()
            detail = info.executer.error_detail
            logging.info("■ The remote process running the task failed:")
            if detail:
                logging.info("\t· Message: %s", detail)
            else:
                logging.info("\t· No error message available.")
        else:
            logging.info("■ %s", description)

    def _update_queue_info(self, is_tty: bool, duration: str) -> None:
        """Update the queue information.

        Prints a message to the user with the number of tasks ahead in the queue
        or that the task is about to start.
        This method replaces the previous message in the terminal.
        """
        #Get status to update the self._tasks_ahead
        self.get_status()
        sys.stdout.write("\r\033[2K")
        sys.stdout.write(
            self._setup_queue_message(is_tty=is_tty, duration=duration))
        sys.stdout.flush()

    def _handle_terminal_status(self, download_std_on_completion: bool,
                                status: models.TaskStatusCode) -> None:
        """Handle a terminal status.

        Downloads the standard files (stdout and stderr) if the task is in a
        terminal status and if the user requested it.
        """
        sys.stdout.flush()
        sys.stdout.write("\r\033[2K")

        if download_std_on_completion:
            self._called_from_wait = True
            out_dir = self.download_outputs(
                filenames=self.STANDARD_OUTPUT_FILES)
            if status == models.TaskStatusCode.FAILED:
                self._print_failed_message(out_dir)

    def wait(self,
             polling_period: int = 1,
             silent_mode: bool = False,
             download_std_on_completion: bool = True) -> models.TaskStatusCode:
        """Wait for the task to complete.

        This method issues requests to the API.

        Args:
            polling_period: How often to poll the API for the task status.
            silent_mode: If True, do not print the task logs (stdout and stderr)
                to the console.
            download_std_on_completion: Request immediate download of the
                standard files (stdout and stderr) after the task completes.

        Returns:
            The final status of the task.
        """
        # TODO: refactor method to make it cleaner
        prev_status = None
        is_tty = sys.stdout.isatty()

        logging.info(
            "Waiting for task %s to complete...\n"
            "Go to https://console.inductiva.ai/tasks/%s for more details.",
            self.id, self.id)

        requires_newline = False
        previous_duration_l = 0

        while True:
            # status = self.get_status()
            task_info = self.get_info()
            status = models.TaskStatusCode(
                task_info.status_history[-1]["status"])
            status_start_time = datetime.datetime.fromisoformat(
                task_info.status_history[-1]["timestamp"])
            description = task_info.status_history[-1].get("description", "")

            now_time = datetime.datetime.now(datetime.timezone.utc)
            duration_timedelta = now_time - status_start_time
            duration_timedelta = max(duration_timedelta, datetime.timedelta(0))
            duration = format_utils.short_timedelta_formatter(
                duration_timedelta)

            if status != prev_status:
                if requires_newline:
                    requires_newline = False
                    sys.stdout.write("\n")
                self._handle_status_change(status, description)

                if (status == models.TaskStatusCode.COMPUTATIONSTARTED) and (
                        not silent_mode):
                    try:
                        self.tail_files(["stdout.txt", "stderr.txt"], 50, True,
                                        sys.stdout)
                    # pylint: disable=broad-except
                    except Exception as _:
                        # Ignore errors while tailing files
                        pass

            # Print timer
            elif (status != models.TaskStatusCode.SUBMITTED and
                  not task_info.is_terminal):

                #clear previous line
                print(" " * previous_duration_l, end="\r")

                duration = f"Duration: {duration}"
                print(duration, end="\r")

                previous_duration_l = len(duration)

            prev_status = status

            #Used to print queue information
            if (status == models.TaskStatusCode.SUBMITTED and
                    self._tasks_ahead is not None):
                requires_newline = True
                self._update_queue_info(is_tty=is_tty, duration=duration)
            #use is_terminal instead of the method to avoid an api call
            #that can make the task status inconsistent
            if self.info.is_terminal:
                self._handle_terminal_status(
                    download_std_on_completion=download_std_on_completion,
                    status=status)
                return status

            time.sleep(polling_period)

    def _validate_task_computation_started(self) -> Tuple[bool, Optional[str]]:
        info = self.get_info()
        if info.is_terminal:
            print(
                f"Task {self.id} has terminated.\n"
                "Access its output using:\n\n"
                f"  inductiva tasks download --id {self.id}",
                file=sys.stderr)
            return False
        if not info.status == "computation-started":
            print(
                f"Task {self.id} has not started yet.\n"
                "Wait for computation to start.",
                file=sys.stderr)
            return False

        return True

    def tail_files(self,
                   tail_files: List[str],
                   lines: int,
                   follow: bool,
                   fout: TextIO = sys.stdout):
        """
        Prints the result of tailing a list of files.

        Args:
            tail_files: A list of files to tail.
            lines: The number of lines to print.
            follow: Whether to keep tailing a file or not. If True, tail_files
                will keep printing the new lines in the selected files as they
                are changed in real time. If False, it will print the tail and
                end.
            fout: The file object to print the result to. Default is stdout.
        """
        return self._run_multiple_streaming_commands([
            lambda filename=filename: self._run_tail_on_machine(
                filename, lines, follow) for filename in tail_files
        ],
                                                     fout=fout)

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
             wait_timeout: Optional[Union[float, int]] = 1,
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
                raise ValueError("Wait timeout must be a positive number"
                                 " or None.")

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
        return self.info.simulator

    def get_storage_path(self) -> str:
        return self.info.storage_path

    def get_output_info(self) -> output_info.TaskOutputInfo:
        """Get information about the output files of the task.

        Returns:
            An instance of the OutputInfo class, which can be used to
            access info about the output archive (number of files, total
            compressed size, total uncompressed size) and information about
            each file (name, size, compressed size). It can also be used to
            print that information in a formatted way.
        """
        archive_info = storage.get_zip_contents(
            path=self.info.storage_output_path, zip_relative_path="artifacts/")

        output_files = [
            output_info.FileInfo(
                name=file_info.name,
                size=file_info.size,
                compressed_size=file_info.compressed_size,
            ) for file_info in archive_info.files
        ]

        return output_info.TaskOutputInfo(
            total_size_bytes=archive_info.size,
            files=output_files,
        )

    def _contains_only_std_files(self, output_dir: pathlib.Path) -> bool:
        """Check if the output archive contains only stdout and stderr files.

        Returns:
            True if the output archive contains only stdout and stderr files,
            False otherwise.
        """
        output_files = list(output_dir.iterdir())
        return all(
            file.name in self.STANDARD_OUTPUT_FILES for file in output_files)

    def _request_download_output_url(self) -> Optional[str]:
        try:
            url = storage.get_signed_urls(paths=[self.info.storage_output_path],
                                          operation="download")[0]
        except exceptions.ApiException as e:
            if not self._called_from_wait:

                if self._status == models.TaskStatusCode.EXECUTERFAILED:
                    logging.info("The remote process running the task failed:")
                    self.get_info()
                    detail = self.info.executer.error_detail
                    if detail:
                        logging.info(" > Message: %s", detail)
                    else:
                        logging.info(" > No error message available.")
                else:
                    # Raise the exception to be handled by the exception handler
                    raise e
            return None
        finally:
            # Reset internal state
            self._called_from_wait = False

        return url

    def _request_download_input_url(self) -> str:
        return storage.get_signed_urls(paths=[self.info.storage_input_path],
                                       operation="download")[0]

    def get_output_url(self) -> Optional[str]:
        """Get a public URL to download the output files of the task.

        Returns:
            The URL to download the output files of the task, or None
            if the
        """
        download_url = self._request_download_output_url()
        if not download_url:
            return None

        logging.info("■ Use the following URL to download the output "
                     "files of you simulation:")
        logging.info(" > %s", download_url)

        return download_url

    def get_input_url(self) -> Optional[str]:
        """Get a public URL to download the input files of the task.

        Returns:
            The URL to download the input files of the task, or None
        """
        download_url = self._request_download_input_url()
        if download_url is None:
            raise RuntimeError(
                "The API did not return a download URL for the task inputs.")

        logging.info("■ Use the following URL to download the input "
                     "files of you simulation:")
        logging.info(" > %s", download_url)

        return download_url

    def _download(
        self,
        filenames: Optional[List[str]],
        dest_dir: Optional[str],
        sub_dir: str,
        uncompress: bool,
        rm_downloaded_zip_archive: bool,
        rm_remote_files: bool,
        zip_name: str,
        request_download_url: Callable,
        download_partial_files: Callable,
    ) -> pathlib.Path:
        self._status = self.get_status()

        download_url = request_download_url()
        if not download_url:
            return None

        logging.debug("\nDownload URL: %s\n", download_url)

        if dest_dir is None:
            dest_dir = self.id
        dest_dir += f"/{sub_dir}/"

        dir_path = files.resolve_output_path(dest_dir)

        if (dir_path.exists() and not self._contains_only_std_files(dir_path)):
            warnings.warn("Path already exists, files may be overwritten.")
        dir_path.mkdir(parents=True, exist_ok=True)

        download_message = "Downloading simulation files to %s..."

        if filenames is self.STANDARD_OUTPUT_FILES:
            download_message = "Downloading stdout and stderr files to %s..."

        if filenames:
            logging.info(download_message, dir_path)
            download_partial_files(download_url, filenames, dir_path)

            # If the user requested a partial download, the full download
            # will be skipped.

            logging.info("Partial download completed to %s.", dir_path)
            return dir_path

        zip_path = dir_path.joinpath(zip_name)
        logging.info(download_message, zip_path)

        pool_manager: urllib3.PoolManager = (
            self._api.api_client.rest_client.pool_manager)
        response = pool_manager.request(
            "GET",
            download_url,
            preload_content=False,
        )

        # use raw urllib3 response instead of the generated client response, to
        # implement our own download logic (with progress bar, first checking
        # the size of the file, etc.)
        data.download_file(response, zip_path)

        if uncompress:
            logging.info("Uncompressing the files to %s...", dir_path)
            data.decompress_zip(zip_path, dir_path)
            if rm_downloaded_zip_archive:
                zip_path.unlink()

        if rm_remote_files:
            self.remove_remote_files()

        if self._status == models.TaskStatusCode.FAILED:
            logging.error(
                "Task %s failed.\n"
                "Please inspect the stdout.txt and stderr.txt files at %s\n"
                "For more information.", self.id, dir_path)

        return dir_path

    def download_outputs(
        self,
        filenames: Optional[List[str]] = None,
        output_dir: Optional[str] = None,
        uncompress: bool = True,
        rm_downloaded_zip_archive: bool = True,
        rm_remote_files: bool = False,
    ) -> Optional[pathlib.Path]:
        """Download output files of the task.

        Args:
            filenames: List of filenames to download. If None or empty, all
                files are downloaded.
            output_dir: Directory where to download the files. If None, the
                files are downloaded to the default directory. The default is
                {inductiva.get_output_dir()}/[{output_dir}|{task_id}]/outputs/.
            uncompress: Whether to uncompress the archive after downloading it.
            rm_downloaded_zip_archive: Whether to remove the archive after
                uncompressing it. If uncompress is False, this argument is
                ignored.
            rm_remote_files: Whether to remove all task files from remote
                storage after the download is complete. Only used if filenames
                is None or empty (i.e., all output files are downloaded).
        """
        return self._download(
            filenames=filenames,
            dest_dir=output_dir,
            sub_dir="outputs",
            uncompress=uncompress,
            rm_downloaded_zip_archive=rm_downloaded_zip_archive,
            rm_remote_files=rm_remote_files,
            zip_name="output.zip",
            request_download_url=self._request_download_output_url,
            download_partial_files=data.download_partial_outputs,
        )

    def download_inputs(
        self,
        filenames: Optional[List[str]] = None,
        input_dir: Optional[str] = None,
        uncompress: bool = True,
        rm_downloaded_zip_archive: bool = True,
        rm_remote_files: bool = False,
    ) -> Optional[pathlib.Path]:
        """Download input files of the task.

        Args:
            filenames: List of filenames to download. If None or empty, all
                files are downloaded.
            input_dir: Directory where to download the files. If None, the
                files are downloaded to the default directory. The default is
                {inductiva.get_output_dir()}/[{input_dir}|{task_id}]/inputs/.
            uncompress: Whether to uncompress the archive after downloading it.
            rm_downloaded_zip_archive: Whether to remove the archive after
                uncompressing it. If uncompress is False, this argument is
                ignored.
            rm_remote_files: Whether to remove all task files from remote
                storage after the download is complete. Only used if filenames
                is None or empty (i.e., all input files are downloaded).
        """
        return self._download(
            filenames=filenames,
            dest_dir=input_dir,
            sub_dir="inputs",
            uncompress=uncompress,
            rm_downloaded_zip_archive=rm_downloaded_zip_archive,
            rm_remote_files=rm_remote_files,
            zip_name="input.zip",
            request_download_url=self._request_download_input_url,
            download_partial_files=data.download_partial_inputs,
        )

    async def _file_operation(self,
                              operation: Operations,
                              formatter: Callable,
                              follow: bool = False,
                              **kwargs):
        """Perform file operations on the task that is currently running.

        Args:
            operation: The operation to perform on the task files.
            **kwargs: Additional arguments for the operation.

        Returns:
            The result of the operation.
        """
        pc = self.file_tracker.create_peer_connection()
        message_queue, end_event = await self.file_tracker.setup_channel(
            pc, operation, follow=follow, **kwargs)
        if not await self.file_tracker.connect_to_task(self._api, pc, self.id):
            yield "Failed to connect to the task."
            return
        while not end_event.is_set():
            message = await message_queue.get()

            if message is None:
                return
            elif message["status"] != "success":
                await self.file_tracker.cleanup()
                yield message["message"]
                return

            yield formatter(message["message"])

        await pc.close()

    async def close_stream(self):
        """Close the stream to the task."""
        if self.file_tracker is not None:
            await self.file_tracker.cleanup()

    def list_files(self) -> Tuple[Optional[str], int]:
        """List the files in the task's working directory.
        
        This method will list the files, in real time, in the task's working
        directory. It will also print the files in a tree-like structure.

        returns:
            A string with the formatted directory listing.
            The return code for the command. 0 if successful, 1 if failed.
        """

        result, return_code = self._run_streaming_command(
            lambda: self._file_operation(
                Operations.LIST, formatter=self._format_directory_listing))

        return result, return_code

    async def _gather_and_consume(self, generators: List[AsyncGenerator],
                                  fout: TextIO):
        """
        Helper method to gather and consume multiple asynchronous generators.
        """
        tasks = [
            asyncio.create_task(self._consume(generator, fout))
            for generator in generators
        ]
        try:
            await asyncio.gather(*tasks)
        except asyncio.CancelledError:
            for task in tasks:
                task.cancel()
            await self.close_stream()

    async def _consume_modified_file(self, generator: AsyncGenerator,
                                     fout: TextIO):
        """
        Consume and write the formatted output from an asynchronous generator to
        a file-like object.

        This function iterates over the provided asynchronous generator, writing 
        each line of output to the specified file-like object.

        Example:
            Most Recent File: /workdir/io9od5da6xh131inmsno0fapm/stdin.txt
            Modification Time: 2025-04-01 09:28:33
            Current Time on Machine: 2025-04-01 09:29:17

            Time Since Last Modification: 0:00:43
        """
        try:
            async for generator_data in generator:

                # Convert timestamps to readable datetime
                most_recent_time = datetime.datetime.fromtimestamp(
                    generator_data["most_recent_timestamp"]).strftime(
                        "%Y-%m-%d %H:%M:%S")
                now_time = datetime.datetime.fromtimestamp(
                    generator_data["now_timestamp"]).strftime(
                        "%Y-%m-%d %H:%M:%S")

                # Print the information
                recent_file = generator_data["most_recent_file"]
                formatted_seconds = format_utils.seconds_formatter(
                    generator_data["time_since_last_mod"])
                print(
                    "\n"
                    f"Most Recent File: {recent_file}\n"
                    f"Modification Time: {most_recent_time}\n"
                    f"Current Time on Machine: {now_time}\n"
                    "\n"
                    f"Time Since Last Modification: {formatted_seconds}",
                    file=fout)
        except asyncio.CancelledError:
            pass

    def _last_modified_file_formatter(self, generator_data: dict) -> str:
        """
        Formats the outputs of the last_modified_file command.
        Args:
            generator_data: The data returned by the last_modified_file
                command.
        """
        # Convert timestamps to readable datetime
        most_recent_time = datetime.datetime.fromtimestamp(
            generator_data["most_recent_timestamp"]).strftime(
                "%Y-%m-%d %H:%M:%S")
        now_time = datetime.datetime.fromtimestamp(
            generator_data["now_timestamp"]).strftime("%Y-%m-%d %H:%M:%S")

        # Print the information
        recent_file = generator_data["most_recent_file"]
        formatted_seconds = format_utils.seconds_formatter(
            generator_data["time_since_last_mod"])
        return ("\n"
                f"Most Recent File: {recent_file}\n"
                f"Modification Time: {most_recent_time}\n"
                f"Current Time on Machine: {now_time}\n"
                "\n"
                f"Time Since Last Modification: {formatted_seconds}")

    def last_modified_file(self):
        """
        Display the last modified file for a given task.

        This function retrieves and prints information about the most recently 
        modified file associated with a specified task. It validates that the 
        task computation has started before proceeding. If the task is invalid 
        or not started, an error message is printed to `stderr`.
        """

        result, return_code = self._run_streaming_command(
            lambda: self._file_operation(
                Operations.LAST_MODIFIED_FILE,
                formatter=self._last_modified_file_formatter,
            ))

        return result, return_code

    async def _run_tail_on_machine(self,
                                   filename: str,
                                   n_lines: int = 10,
                                   follow=False):
        """Get the last n_lines lines of a 
        file in the task's working directory."""

        def formatter(message):
            return self._format_list_of_lines(message,
                                              filename,
                                              endl="\n",
                                              header=not follow)

        async for lines in self._file_operation(Operations.TAIL,
                                                formatter=formatter,
                                                filename=filename,
                                                lines=n_lines,
                                                follow=follow):
            yield lines

    async def _consume(self, generator: AsyncGenerator, fout: TextIO):
        """
        Consume and write the output from an asynchronous generator to a
        file-like object.

        This function iterates over the provided asynchronous generator, writing 
        each line of output to the specified file-like object.
        """
        try:
            async for lines in generator:
                print(lines, file=fout, end="", flush=True)
        except asyncio.CancelledError:
            pass

    def _top(self) -> Tuple[Optional[str], int]:
        """Prints the result of the `top -b -H -n 1` command.
    
        This command will list the processes and threads (-H) in batch mode
        (-b).
        This command will run only once (-n 1) instead of running continuously.
        The result is an instant snapshot of the machine CPU and RAM metrics.

        Returns:
            A string with the formatted directory listing. 
            The return code for the command. 0 if successful, 1 if failed.
        """
        result, return_code = self._run_streaming_command(
            lambda: self._file_operation(
                Operations.TOP, formatter=lambda _: _, follow=False))

        return result, return_code

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
        if self.info.is_terminal:
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

    def remove_remote_files(self, verbose: bool = True) -> bool:
        """Removes all files associated with the task from remote storage.

        Returns:
            True if the files were removed successfully, False otherwise.
        """
        if verbose:
            logging.info(
                "Removing files from remote storage for task %s...",
                self.id,
            )
        try:
            # TODO: rename the function to a more generic name
            storage.remove_workspace(remote_dir=self.id)
            if verbose:
                logging.info("Remote task files removed successfully.")
        except exceptions.ApiException as e:
            logging.error("An error occurred while removing the files:")
            logging.error(" > %s", json.loads(e.body)["detail"])
            return False
        return True

    def _get_summary(self) -> str:
        """Get a formatted summary of the task. This method caches the
        summary."""
        info: TaskInfo = self.get_info()

        # Update the duration metrics if the task is still running, otherwise
        # the cached values will be used
        info.time_metrics.total_seconds.value = self.get_total_time(cached=True)
        info.time_metrics.computation_seconds.value = \
            self.get_computation_time(cached=True)

        self._summary = str(info)
        return self._summary

    def _run_multiple_streaming_commands(
            self,
            generator_factories: List[Callable[[], AsyncGenerator]],
            fout: TextIO = sys.stdout):
        if not self._validate_task_computation_started():
            return 1

        if inductiva.is_notebook():
            nest_asyncio.apply()

        asyncio.run(
            self._gather_and_consume([gen() for gen in generator_factories],
                                     fout))
        return 0

    def _run_streaming_command(
        self, generator_factory: Callable[[], AsyncGenerator]
    ) -> Tuple[Optional[str], int]:
        if not self._validate_task_computation_started():
            return None, 1

        if inductiva.is_notebook():
            nest_asyncio.apply()

        buffer = io.StringIO()

        asyncio.run(self._gather_and_consume([generator_factory()], buffer))

        return buffer.getvalue(), 0

    @property
    def summary(self) -> str:
        """It returns cached information about the task summary."""
        if self._summary is None:
            return self._get_summary()
        return self._summary

    def print_summary(self, fhandle=sys.stdout):
        print(self._get_summary(), file=fhandle)

    def set_metadata(self, metadata: Dict[str, str]):
        """Set metadata for the task.
        
        Metadata is stored as key-value pairs, where both
        keys and values must be strings.
        Metadata can be useful for categorizing, searching,
        and filtering tasks.

        Example usage:
            task = simulator.run(...)
            # Add experiment information to the task
            task.set_metadata({
                "study": "study_1",
                "experiment": "experiment_1",
                "description": "This is a test experiment",
                "parameters": "param1=1,param2=2",
            })

        Args:
            metadata: A dictionary with the metadata to set.
        """
        # Validate metadata
        if not isinstance(metadata, Dict):
            raise TypeError("Metadata must be a dictionary.")

        for key, value in metadata.items():
            if not isinstance(key, str) or not isinstance(value, str):
                raise TypeError("Metadata keys and values must be strings.")
            if not key or not value:
                raise ValueError(
                    "Metadata keys and values cannot be empty strings.")

        self._api.set_metadata(path_params={"task_id": self.id}, body=metadata)

    def get_metadata(self) -> Dict[str, str]:
        """Get the metadata associated with the task.
            
        Returns:
            A dictionary with the custom metadata previously set on this task.
        """
        response = self._api.get_metadata(path_params={"task_id": self.id})

        return dict(response.body)
