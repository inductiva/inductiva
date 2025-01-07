"""Methods that interact with the lower-level inductiva-web-api-client.

The relevant function for usage outside of this file is invoke_async_api().
Check the demos directory for examples on how it is used.
"""
import os
import time
import tqdm
import tqdm.utils
import signal
import urllib3
import decimal
from contextlib import contextmanager
from typing import Any, Dict, List, Optional, Type

import logging

import inductiva
from inductiva.client import ApiClient, ApiException, Configuration
from inductiva.client.apis.tags.tasks_api import TasksApi
from inductiva.client.models import TaskRequest, TaskStatus, TaskSubmittedInfo
from inductiva import types, constants
from inductiva.utils.data import (get_validate_request_params, pack_input)
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


def prepare_input(task_id, original_params, type_annotations):
    sim_dir = original_params["sim_dir"]
    # If the input directory is empty, do not zip it
    # still need to zip the input parameters though
    if sim_dir:
        inputs_size = files.get_path_size(sim_dir)
        logging.info("Preparing upload of the local input directory %s (%s).",
                     sim_dir, format_utils.bytes_formatter(inputs_size))

        if os.path.isfile(os.path.join(sim_dir, constants.TASK_OUTPUT_ZIP)):
            raise ValueError(
                f"Invalid file name: '{constants.TASK_OUTPUT_ZIP}'")

    input_zip_path = pack_input(
        params=original_params,
        type_annotations=type_annotations,
        zip_name=task_id,
    )

    zip_file_size = os.path.getsize(input_zip_path)
    logging.info("Input archive size: %s",
                 format_utils.bytes_formatter(zip_file_size))

    logging.info("Uploading input archive...")

    return input_zip_path, zip_file_size


def get_upload_url(
    api_endpoint,
    query_params: Dict[str, str] = None,
    path_params: Dict[str, str] = None,
):
    """
    Fetches the upload URL using the specified method.
    """
    params = {}
    if query_params is not None:
        params["query_params"] = query_params
    if path_params is not None:
        params["path_params"] = path_params

    return api_endpoint(**params).body


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


def notify_upload_complete(api_endpoint,
                           query_params: Dict[str, str] = None,
                           path_params: Dict[str, str] = None):
    """
    Notifies the API that the file upload is complete.
    """
    params = {}
    if query_params is not None:
        params["query_params"] = query_params
    if path_params is not None:
        params["path_params"] = path_params

    api_endpoint(**params)


def upload_input(api_instance: TasksApi, task_id, original_params,
                 type_annotations):
    """Uploads the inputs of a given task to the API.

    Args:
        api_instance: Instance of TasksApi used to send necessary requests.
        task_id: ID of the task.
        original_params: Params of the request passed by the user.
        type_annotations: Annotations of the params' types.
    """

    input_zip_path, zip_file_size = prepare_input(task_id, original_params,
                                                  type_annotations)

    api_response = get_upload_url(
        api_instance.get_input_upload_url,
        path_params={"task_id": task_id},
    )

    with tqdm.tqdm(total=zip_file_size,
                   unit="B",
                   unit_scale=True,
                   unit_divisor=1000) as progress_bar:

        method = api_response["method"]
        url = api_response["url"]

        upload_file(api_instance, input_zip_path, method, url, progress_bar)

        notify_upload_complete(
            api_instance.notify_input_uploaded,
            path_params={"task_id": task_id},
        )

    logging.info("Local input directory successfully uploaded.")
    logging.info("")
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


def task_info_str(
    task_id,
    params,
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

    local_input_dir = params["sim_dir"]
    info_str += (f"\t· Local input directory: {local_input_dir}\n"
                 "\t· Submitting to the following computational resources:\n")
    if resource_pool is not None:
        info_str += f" \t\t· {resource_pool}\n"
    else:
        machine_type = constants.DEFAULT_QUEUE_MACHINE_TYPE
        info_str += f" \t\t· Default queue with {machine_type} machines.\n"
        ttl_seconds = task_submitted_info.get("time_to_live_seconds")
        if ttl_seconds is not None and isinstance(ttl_seconds, decimal.Decimal):
            ttl_seconds = format_utils.seconds_formatter(ttl_seconds)
            info_str += (f" \t\t· Task will be killed after the computation "
                         f"time exceeds {ttl_seconds} (h:m:s).\n")
    info_str += "\n"
    return info_str


def submit_task(api_instance,
                simulator,
                request_params,
                resource_pool,
                storage_path_prefix,
                params,
                type_annotations,
                resubmit_on_preemption: bool = False,
                container_image: Optional[str] = None,
                simulator_name_alias: Optional[str] = None,
                simulator_obj=None,
                input_resources: Optional[List[str]] = None):
    """Submit a task and send input files to the API."""
    resource_pool_id = resource_pool.id

    current_project = inductiva.projects.get_current_project()
    if current_project is not None:
        if not current_project.opened:
            raise RuntimeError("Trying to submit a task to a closed project.")
        current_project = current_project.name

    if not input_resources:
        input_resources = []
    task_request = TaskRequest(simulator=simulator,
                               params=request_params,
                               project=current_project,
                               resource_pool=resource_pool_id,
                               container_image=container_image,
                               storage_path_prefix=storage_path_prefix,
                               simulator_name_alias=simulator_name_alias,
                               resubmit_on_preemption=resubmit_on_preemption,
                               input_resources=input_resources)

    task_submitted_info = submit_request(
        api_instance=api_instance,
        request=task_request,
    )

    task_id = task_submitted_info["id"]
    logging.info(
        task_info_str(
            task_id,
            params,
            resource_pool,
            simulator_obj,
            task_submitted_info,
        ))

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
                     simulator_obj=None,
                     simulator_name_alias=None,
                     input_resources: Optional[List[str]] = None) -> str:
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
        input_resources: Additional input files that will be copied to the
                simulation from a bucket or from another task output.
    Return:
        Returns the task id.
    """

    request_params = get_validate_request_params(
        original_params=params,
        type_annotations=type_annotations,
    )

    with get_client() as client:
        api_instance = TasksApi(client)

        task_id = submit_task(api_instance=api_instance,
                              simulator=simulator,
                              simulator_name_alias=simulator_name_alias,
                              request_params=request_params,
                              resource_pool=resource_pool,
                              storage_path_prefix=storage_path_prefix,
                              params=params,
                              container_image=container_image,
                              type_annotations=type_annotations,
                              resubmit_on_preemption=resubmit_on_preemption,
                              simulator_obj=simulator_obj,
                              input_resources=input_resources)

    return task_id
