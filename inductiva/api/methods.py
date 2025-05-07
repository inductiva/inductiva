"""Methods that interact with the lower-level inductiva-web-api-client.

The relevant function for usage outside of this file is submit_task().
Check the demos directory for examples on how it is used.
"""
import os
import sys
import time
import tqdm
import tqdm.utils
import signal
import urllib3
import decimal
from contextlib import contextmanager
from typing import List, Optional

import logging

import inductiva
from inductiva.client import ApiClient, ApiException, Configuration
from inductiva.client.apis.tags.tasks_api import TasksApi
from inductiva.client.models import (TaskRequest, TaskStatus, TaskSubmittedInfo,
                                     CompressionMethod)
from inductiva import constants, storage
from inductiva.utils import format_utils, files


def get_api_config() -> Configuration:
    """Returns an API configuration object."""
    api_key = inductiva.get_api_key()

    api_config = Configuration(host=inductiva.api_url)
    api_config.api_key["APIKeyHeader"] = api_key

    return api_config


def get_client(api_config: Optional[Configuration] = None) -> ApiClient:
    """Returns an ApiClient instance."""

    if api_config is None:
        api_config = get_api_config()

    client = ApiClient(api_config)

    client.user_agent = inductiva.get_api_agent()

    return client


def submit_request(task_api_instance: TasksApi,
                   request: TaskRequest) -> TaskSubmittedInfo:
    """Submits a task request to the API.

    Args:
        api_instance: Instance of TasksApi used to send necessary requests.
    Return:
        Returns the body of the HTTP response.
        Contains two fields, "id" and "status".
    """
    # Submit task using provided or temporary instance
    api_response = task_api_instance.submit_task(body=request)
    logging.debug("Request status: %s", api_response.body["status"])

    return api_response.body


def prepare_input(task_id, input_dir, params):
    """Prepare the input files for a task submission."""

    # If the input directory is empty, do not zip it
    # still need to zip the input params parameters though
    if input_dir:
        inputs_size = files.get_path_size(input_dir)
        logging.info("Preparing upload of the local input directory %s (%s).",
                     input_dir, format_utils.bytes_formatter(inputs_size))

        if os.path.isfile(os.path.join(input_dir, constants.TASK_OUTPUT_ZIP)):
            raise ValueError(
                f"Invalid file name: '{constants.TASK_OUTPUT_ZIP}'")

    input_zip_path = inductiva.utils.data.pack_input(
        input_dir,
        params,
        zip_name=task_id,
    )

    zip_file_size = os.path.getsize(input_zip_path)
    logging.info("Input archive size: %s",
                 format_utils.bytes_formatter(zip_file_size))

    logging.info("Uploading input archive...")

    return input_zip_path, zip_file_size


def upload_file(api_instance: ApiClient, input_path: str, method: str, url: str,
                progress_bar: tqdm):
    """
    Handles the upload of a file, updating the provided progress bar.
    """
    headers = {"Content-Type": "application/octet-stream"}

    with open(input_path, "rb") as zip_fp:
        wrapped_file = tqdm.utils.CallbackIOWrapper(progress_bar.update, zip_fp,
                                                    "read")
        pool_manager: urllib3.PoolManager = (
            api_instance.api_client.rest_client.pool_manager)
        resp = pool_manager.request(method,
                                    url,
                                    body=wrapped_file,
                                    headers=headers)
        if resp.status != 200:
            raise ApiException(status=resp.status, reason=resp.reason)


def upload_input(api_instance: TasksApi, input_dir, params, task_id,
                 storage_path_prefix):
    """Uploads the inputs of a given task to the API.

    Args:
        api_instance: Instance of TasksApi used to send necessary requests.
        task_id: ID of the task.
        input_dir: Directory containing the input files to be uploaded.
        params: Additional parameters to be sent to the API.
        storage_path_prefix: Path to the storage bucket.
        """
    input_zip_path = None

    try:
        input_zip_path, zip_file_size = prepare_input(task_id, input_dir,
                                                      params)

        remote_input_zip_path = f"{storage_path_prefix}/{task_id}/input.zip"
        url = storage.get_signed_urls(
            paths=[remote_input_zip_path],
            operation="upload",
        )[0]

        with tqdm.tqdm(total=zip_file_size,
                       unit="B",
                       unit_scale=True,
                       unit_divisor=1000) as progress_bar:
            upload_file(api_instance, input_zip_path, "PUT", url, progress_bar)
            api_instance.notify_input_uploaded(path_params={"task_id": task_id})
        logging.info("Local input directory successfully uploaded.")
        logging.info("")

    finally:
        if input_zip_path:
            os.remove(input_zip_path)


def block_until_finish(api_instance: TasksApi, task_id: str) -> str:
    """Block until a task executing remotely finishes execution.

    Args:
        api_instance: Instance of TasksApi used to send necessary requests.
        task_id: ID of the task to wait for.

    Returns:
        Returns info related to the task, containing two fields,
        "id" and "status".
    """
    logging.debug("Blocking until task is finished ...")
    return block_until_status_is(api_instance, task_id, {"success", "failed"})


def kill_task(api_instance: TasksApi, task_id: str):
    """Kill a task that is executing remotely.

    The function sends a kill request to the API.

    Args:
        api_instance: Instance of TasksApi used to send necessary requests.
        task_id: ID of the task to kill.
   """
    logging.debug("Sending kill task request ...")
    api_instance.kill_task(path_params={"task_id": task_id},)
    logging.info("Task with ID %s was terminated.", task_id)


def get_task_status(api_instance: TasksApi, task_id: str) -> TaskStatus:
    """Check the status of a task."""

    api_response = api_instance.get_task_status(
        path_params={"task_id": task_id})

    status = api_response.body["status"]

    return status


def block_until_status_is(api_instance: TasksApi,
                          task_id,
                          desired_status,
                          sleep_secs=0.5):
    """Block until the status of a task becomes the desired status.

    Args:
        api_instance: Instance of TasksApi used to send necessary requests.
        task_id: ID of the task to wait for.
        desired_status: Task status to wait for.
        sleep_secs: Polling interval.

    Returns:
        Returns info related to the task, containing two fields,
    """
    prev_status = None

    while True:
        try:
            status = get_task_status(api_instance, task_id)
            logging.debug("Task status is %s", status)
        except ApiException as e:
            raise e

        if status != prev_status:
            if status == "submitted":
                logging.info("Waiting for resources...")

            prev_status = status

        # If status reaches the desired status, then stop polling
        if status in desired_status:
            break

        time.sleep(sleep_secs)

    return status


def _configure_sigint_handler(handler):
    if not handler:
        return None

    try:
        return signal.signal(signal.SIGINT, handler)
    except ValueError as e:
        logging.warning("Custom SIGINT handler not configured: %s", e)
        # If the signal is not supported, ignore it
        pass


@contextmanager
def blocking_task_context(api_instance: TasksApi,
                          task_id: str,
                          action_str: str = "action"):
    """Context to handle execution of a blocking task.

    The context handles exceptions and the SIGINT signal, issuing a request
    to the API to kill the executing task.
    Info on the implementation:
    - https://docs.python.org/3/library/contextlib.html

    Args:
        api_instance: Instance of TasksApi used to send necessary requests.
        task_id: ID of the task being executed.
    """
    # Other imported modules can make changes to the SIGINT handler, so we set
    # the default int handler to make sure that KeyboardInterrupt is raised
    # if the user presses Ctrl+C.
    original_sig = _configure_sigint_handler(signal.default_int_handler)

    try:
        yield None
    except Exception as err:
        logging.info("Caught exception: terminating blocking task...")
        kill_task(api_instance, task_id)
        raise err
    except KeyboardInterrupt:
        logging.info("Caught SIGINT: %s interrupted by user.", action_str)
        kill_task(api_instance, task_id)
        sys.exit(1)
    finally:
        # Reset original SIGINT handler
        _configure_sigint_handler(original_sig)


def task_info_str(
    task_id,
    local_input_dir,
    resource_pool,
    simulator,
    task_submitted_info: TaskSubmittedInfo,
) -> str:
    """Generate a string with the main components of a task submission."""

    info_str = ("■ Task Information:\n"
                f"\t· ID:                    {task_id}\n")
    if simulator is not None:
        info_str += (f"\t· Simulator:             {simulator.name}\n"
                     f"\t· Version:               {simulator.version}\n"
                     f"\t· Image:                 {simulator.image_uri}\n")

    info_str += (f"\t· Local input directory: {local_input_dir}\n"
                 "\t· Submitting to the following computational resources:\n")
    info_str += f" \t\t· {resource_pool}\n"

    if task_submitted_info is not None:
        ttl_seconds = task_submitted_info.get("time_to_live_seconds")
        if ttl_seconds is not None and isinstance(ttl_seconds, decimal.Decimal):
            ttl_seconds = format_utils.seconds_formatter(ttl_seconds)
            info_str += (f" \t\t· Task will be killed after the computation "
                         f"time exceeds {ttl_seconds} (h:m:s).\n")
    info_str += "\n"
    return info_str


def submit_task(simulator,
                input_dir,
                machine_group,
                params,
                storage_path_prefix,
                resubmit_on_preemption: bool = False,
                container_image: Optional[str] = None,
                simulator_name_alias: Optional[str] = None,
                simulator_obj=None,
                remote_assets: Optional[List[str]] = None,
                project_name: Optional[str] = None):
    """Submit a task and send input files to the API.

    Args:
        simulator: The simulator to use
        input_dir: Directory containing the input files to be uploaded.
        machine_group: Group of machines with a queue to submit the task to.
        params: Additional parameters to pass to the simulator.
        storage_path_prefix: Path prefix for storing simulation data
        resubmit_on_preemption (bool): Resubmit task for execution when
                previous execution attempts were preempted. Only applicable when
                using a preemptible resource, i.e., resource instantiated with
                `spot=True`.
        container_image: The container image to use for the simulation
            Example: container_image="docker://inductiva/kutu:xbeach_v1.23_dev"
        simulator_name_alias: Optional alias name for the simulator
        simulator_obj: Optional simulator object with additional configuration
        remote_assets: Additional input files that will be copied to the
                simulation from a bucket or from another task output.
        project: Name of the project to which the task will be
                assigned. If None, the task will be assigned to
                the default project.
    Return:
        Returns the task id.
    """

    if not remote_assets:
        remote_assets = []

    stream_zip = params.pop("stream_zip", True)
    compress_with = params.pop("compress_with", CompressionMethod.SEVEN_Z)

    task_request = TaskRequest(simulator=simulator,
                               params=params,
                               project=project_name,
                               resource_pool=machine_group.id,
                               container_image=container_image,
                               storage_path_prefix=storage_path_prefix,
                               simulator_name_alias=simulator_name_alias,
                               resubmit_on_preemption=resubmit_on_preemption,
                               input_resources=remote_assets,
                               stream_zip=stream_zip,
                               compress_with=compress_with)

    # Create an instance of the TasksApi class
    task_api_instance = TasksApi(get_client())

    # Submit task via the "POST task/submit" endpoint.
    # HTTP status code 400 informs the requested method is invalid.
    # HTTP status code 403 informs that the user is not authorized.
    task_submitted_info = submit_request(
        task_api_instance=task_api_instance,
        request=task_request,
    )

    task_id = task_submitted_info["id"]
    logging.info(
        task_info_str(
            task_id,
            input_dir,
            machine_group,
            simulator_obj,
            task_submitted_info,
        ))

    # If the status returned by the previous HTTP request is "pending-input",
    #  ZIP inputs and send them via "POST task/{task_id}/input".
    if task_submitted_info["status"] == "pending-input":
        # Use the blocking task context
        with blocking_task_context(task_api_instance, task_id, "input upload"):
            upload_input(
                api_instance=task_api_instance,
                input_dir=input_dir,
                params=params,
                task_id=task_id,
                storage_path_prefix=storage_path_prefix,
            )

    # Return task_id and leaves the simulation on the queue until resources
    # become available.
    return task_id
