"""Methods that interact with the lower-level inductiva-web-api-client.

The relevant function for usage outside of this file is invoke_async_api().
Check the demos directory for examples on how it is used.
"""
import os
import time
import tqdm
import tqdm.utils
import signal
import pathlib
import urllib3
import decimal
from contextlib import contextmanager
from typing import Any, Dict, List, Optional, Tuple, Type

import logging

import inductiva
from inductiva.client import ApiClient, ApiException, Configuration
from inductiva.client.apis.tags.tasks_api import TasksApi
from inductiva.client.models import TaskRequest, TaskStatus, TaskSubmittedInfo
from inductiva import types, constants
from inductiva.utils.data import (extract_output, get_validate_request_params,
                                  pack_input)
from inductiva.utils import format_utils, files


def get_api_config() -> Configuration:
    """Returns an API configuration object."""
    api_key = inductiva.get_api_key()

    api_config = Configuration(host=inductiva.api_url)
    api_config.api_key["APIKeyHeader"] = api_key

    return api_config


def get_client() -> ApiClient:
    """Returns an ApiClient instance."""
    api_config = get_api_config()

    return ApiClient(api_config)


def submit_request(api_instance: TasksApi,
                   request: TaskRequest) -> TaskSubmittedInfo:
    """Submits a task request to the API.

    Args:
        api_instance: Instance of TasksApi used to send necessary requests.
    Return:
        Returns the body of the HTTP response.
        Contains two fields, "id" and "status".
    """

    api_response = api_instance.submit_task(body=request)

    logging.debug("Request status: %s", api_response.body["status"])

    return api_response.body


def upload_input(api_instance: TasksApi, task_id, original_params,
                 type_annotations):
    """Uploads the inputs of a given task to the API.

    Args:
        api_instance: Instance of TasksApi used to send necessary requests.
        task_id: ID of the task.
        original_params: Params of the request passed by the user.
        type_annotations: Annotations of the params' types.
    """

    inputs_size = files.get_path_size(original_params["sim_dir"])
    logging.info("Preparing upload of the local input directory %s (%s).",
                 original_params["sim_dir"],
                 format_utils.bytes_formatter(inputs_size))
    input_zip_path = pack_input(
        params=original_params,
        type_annotations=type_annotations,
        zip_name=task_id,
    )

    zip_file_size = os.path.getsize(input_zip_path)
    logging.info("Input archive size: %s",
                 format_utils.bytes_formatter(zip_file_size))

    logging.info("Uploading input archive...")

    try:
        api_response = api_instance.get_input_upload_url(
            path_params={"task_id": task_id})

        method = api_response.body["method"]
        url = api_response.body["url"]
        file_server_available = bool(api_response.body["file_server_available"])

        headers = {"Content-Type": "application/octet-stream"}

        if file_server_available is False:
            headers["X-API-Key"] = (
                api_instance.api_client.configuration.api_key["APIKeyHeader"])

        logging.debug("Upload URL: %s", url)

        with open(input_zip_path, "rb") as zip_fp, tqdm.tqdm(
                total=zip_file_size,
                unit="B",
                unit_scale=True,
                unit_divisor=1000,
        ) as progress_bar:
            # Wrap the file object so that the progress bar is updated every
            # time a chunk is read.
            wrapped_file = tqdm.utils.CallbackIOWrapper(
                progress_bar.update,
                zip_fp,
                "read",
            )

            # Use the pool_manager from the API client to send the request
            # instead of using the generated client. This is because the
            # generated client implementation does not support streaming
            # the file and does not provide a way to update the progress bar.
            pool_manager: urllib3.PoolManager = (
                api_instance.api_client.rest_client.pool_manager)

            resp = pool_manager.request(
                method,
                url,
                body=wrapped_file,
                headers=headers,
            )
            if resp.status != 200:
                raise ApiException(
                    status=resp.status,
                    reason=resp.reason,
                )

            api_response = api_instance.notify_input_uploaded(
                path_params={"task_id": task_id})
    except ApiException as e:
        logging.exception("Exception while uploading local input directory: %s",
                          e)
        raise e

    logging.info("Local input directory successfully uploaded.")
    logging.info("")
    os.remove(input_zip_path)


def download_output(
        api_instance: TasksApi,
        task_id,
        output_dir: Optional[str] = None) -> Tuple[List, pathlib.Path]:
    """Downloads the output of a given task from the API.

    Args:
        api_instance: Instance of TasksApi used to send necessary requests.
        task_id: ID of the task.

    Return:
        Downloads and extracts the data to the client.
    """

    if output_dir is None:
        output_dir = os.path.join(inductiva.get_output_dir(), task_id)

    logging.info("Downloading the task %s outputs to %s...", task_id,
                 output_dir)
    try:
        api_response = api_instance.download_task_output(
            path_params={"task_id": task_id},
            stream=True,
        )
    except ApiException as e:
        raise e

    logging.debug("Downloaded output to %s", api_response.body.name)

    result_list = extract_output(api_response.body.name, output_dir)
    logging.info("Task %s output successfully downloaded to %s.", task_id,
                 output_dir)

    return result_list, pathlib.Path(output_dir)


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
    logging.info("Task terminated.")


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


@contextmanager
def blocking_task_context(api_instance: TasksApi, task_id):
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
    original_sig = signal.signal(signal.SIGINT, signal.default_int_handler)

    try:
        yield None
    except Exception as err:
        logging.info("Caught exception: terminating blocking task...")
        kill_task(api_instance, task_id)
        raise err
    except KeyboardInterrupt as err:
        logging.info("Caught SIGINT: terminating blocking task...")
        kill_task(api_instance, task_id)
        raise err
    finally:
        # Reset original SIGINT handler
        signal.signal(signal.SIGINT, original_sig)


def log_task_info(
    task_id,
    params,
    resource_pool,
    simulator,
    task_submitted_info: TaskSubmittedInfo,
):
    """Logging the main components of a task submission."""

    logging.info("■ Task Information:")
    logging.info("\t· ID:                    %s", task_id)
    if simulator is not None:
        logging.info("\t· Simulator:             %s", simulator.name)
        logging.info("\t· Version:               %s", simulator.version)
        logging.info("\t· Image:                 %s", simulator.image_uri)

    logging.info("\t· Local input directory: %s", params["sim_dir"])
    logging.info("\t· Submitting to the following computational resources:")
    if resource_pool is not None:
        logging.info(" \t\t· %s", resource_pool)
    else:
        logging.info(" \t\t· Default queue with %s machines.",
                     constants.DEFAULT_QUEUE_MACHINE_TYPE)
        ttl_seconds = task_submitted_info.get("time_to_live_seconds")
        if ttl_seconds is not None and isinstance(ttl_seconds, decimal.Decimal):
            logging.info(
                (" \t\t· Task will be killed after the computation time "
                 "exceeds %s (h:m:s)."),
                format_utils.seconds_formatter(ttl_seconds),
            )
    logging.info("")


def submit_task(api_instance,
                simulator,
                request_params,
                resource_pool,
                storage_path_prefix,
                params,
                type_annotations,
                resubmit_on_preemption: bool = False,
                container_image: Optional[str] = None,
                simulator_obj=None):
    """Submit a task and send input files to the API."""
    resource_pool_id = resource_pool.id

    current_project = inductiva.projects.get_current_project()
    if current_project is not None:
        if not current_project.opened:
            raise RuntimeError("Trying to submit a task to a closed project.")
        current_project = current_project.name

    task_request = TaskRequest(simulator=simulator,
                               params=request_params,
                               project=current_project,
                               resource_pool=resource_pool_id,
                               container_image=container_image,
                               storage_path_prefix=storage_path_prefix,
                               resubmit_on_preemption=resubmit_on_preemption)

    task_submitted_info = submit_request(
        api_instance=api_instance,
        request=task_request,
    )

    task_id = task_submitted_info["id"]
    log_task_info(
        task_id,
        params,
        resource_pool,
        simulator_obj,
        task_submitted_info,
    )

    if task_submitted_info["status"] == "pending-input":

        upload_input(
            api_instance=api_instance,
            original_params=params,
            task_id=task_id,
            type_annotations=type_annotations,
        )

    return task_id


def invoke_async_api(simulator: str,
                     params,
                     type_annotations: Dict[Any, Type],
                     resource_pool: types.ComputationalResources,
                     storage_path_prefix: Optional[str] = "",
                     container_image: Optional[str] = None,
                     resubmit_on_preemption: bool = False,
                     simulator_obj=None) -> str:
    """Perform a task asyc and remotely via Inductiva's Web API.

    Submits a simulation async to the API and returns the task id.
    The flow is summarized as follows:
        1. Transform request params into the params used to
        validate permission to execute the request.
        2. Submit task via the "POST task/submit" endpoint.
            Note: HTTP status code 400 informs the requested method is invalid,
                and HTTP status code 403 informs that the user is not authorized
                to post such request.
        3. If the status returned by the previous HTTP request is
            "pending-input", ZIP inputs and send them via
            "POST task/{task_id}/input".
        4. Return task_id and leaves the simulation on the queue until resources
            become available.

    Args:
        request: Request sent to the API for validation.
        input_dir: Directory containing the input files to be uploaded.
        container_image: The container image to use for the simulation
            Example: container_image="docker://inductiva/kutu:xbeach_v1.23_dev"
        resubmit_on_preemption (bool): Resubmit task for execution when
                previous execution attempts were preempted. Only applicable when
                using a preemptible resource, i.e., resource instantiated with
                `spot=True`.
    Return:
        Returns the task id.
    """

    api_config = get_api_config()

    request_params = get_validate_request_params(
        original_params=params,
        type_annotations=type_annotations,
    )

    with ApiClient(api_config) as client:
        api_instance = TasksApi(client)

        task_id = submit_task(api_instance=api_instance,
                              simulator=simulator,
                              request_params=request_params,
                              resource_pool=resource_pool,
                              storage_path_prefix=storage_path_prefix,
                              params=params,
                              container_image=container_image,
                              type_annotations=type_annotations,
                              resubmit_on_preemption=resubmit_on_preemption,
                              simulator_obj=simulator_obj)

    return task_id
