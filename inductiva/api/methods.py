"""Methods that interact with the lower-level inductiva-web-api-client.

The relevant function for usage outside of this file is invoke_async_api().
Check the demos directory for examples on how it is used.
"""
import os
import pathlib
import signal
import time
from urllib3.exceptions import MaxRetryError, NewConnectionError
from contextlib import contextmanager
from typing import Any, Dict, List, Optional, Tuple, Type

from absl import logging

import inductiva
from inductiva.client import ApiClient, ApiException, Configuration
from inductiva.client.apis.tags.tasks_api import TasksApi
from inductiva.client.apis.tags.version_api import VersionApi
from inductiva.client.models import (BodyUploadTaskInput, TaskRequest,
                                     TaskStatus)
from inductiva import types, constants
from inductiva.utils.data import (extract_output, get_validate_request_params,
                                  pack_input)
from inductiva.utils import format_utils


def validate_api_key(api_key: Optional[str]) -> Configuration:
    """Validates the API key and returns API configuration"""
    if inductiva.api_key is None:
        # pylint: disable=line-too-long
        raise ValueError(
            "No API Key specified. "
            "Please set the INDUCTIVA_API_KEY environment variable.\n"
            "More infomation at:"
            "https://inductiva-research-labs-inductiva.readthedocs-hosted.com/en/latest/Getting%20Started.html"
        )
    # pylint: enable=line-too-long

    # Perform version check only on first invocation
    if not hasattr(validate_api_key, "version_checked"):
        compare_client_and_backend_versions(inductiva.__version__)
        validate_api_key.version_checked = True

    api_config = Configuration(host=inductiva.api_url)
    api_config.api_key["APIKeyHeader"] = api_key

    return api_config


def get_client() -> ApiClient:
    """Returns an ApiClient instance."""
    api_config = validate_api_key(inductiva.api_key)

    return ApiClient(api_config)


def submit_request(api_instance: TasksApi, request: TaskRequest) -> TaskStatus:
    """Submits a task request to the API.

    Args:
        api_instance: Instance of TasksApi used to send necessary requests.
    Return:
        Returns the body of the HTTP response.
        Contains two fields, "id" and "status".
    """

    try:
        api_response = api_instance.submit_task(body=request)
    except ApiException as e:
        logging.exception("Exception when calling TasksApi->submit_task: %s", e)
        raise e

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

    inputs_size = os.path.getsize(original_params["sim_dir"])
    logging.info("Preparing upload of the local input directory %s (%s).",
                 original_params["sim_dir"],
                 format_utils.bytes_formatter(inputs_size))
    input_zip_path = pack_input(
        params=original_params,
        type_annotations=type_annotations,
        zip_name=task_id,
    )

    try:
        with open(input_zip_path, "rb") as zip_fp:
            _ = api_instance.upload_task_input(
                path_params={"task_id": task_id},
                body=BodyUploadTaskInput(file=zip_fp),
            )
    except ApiException as e:
        logging.exception(
            "Exception when calling TasksApi->upload_task_inputs: %s", e)
        raise e

    logging.info("Local input directory successfully uploaded.")

    os.remove(input_zip_path)


def download_output(
        api_instance: TasksApi,
        task_id,
        output_dir: Optional[types.Path] = None) -> Tuple[List, pathlib.Path]:
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


def log_task_info(task_id, method_name, params, resource_pool):
    """Logging the main components of a task submission."""

    logging.info("Task Information:")
    logging.info("> ID:                    %s", task_id)
    logging.info("> Method:                %s", method_name.split(".")[1])
    logging.info("> Local input directory: %s", params["sim_dir"])
    logging.info("> Submitting to the following computational resources:")
    if resource_pool is not None:
        logging.info(" >> %s", resource_pool)
    else:
        logging.info(" >> Default queue with %s machines.",
                     constants.DEFAULT_QUEUE_MACHINE_TYPE)


def submit_task(api_instance, method_name, request_params, resource_pool,
                storage_path_prefix, params, type_annotations):
    """Submit a task and send input files to the API."""

    resource_pool_id = None
    if resource_pool is not None:
        resource_pool_id = resource_pool.id

    task_request = TaskRequest(method=method_name,
                               params=request_params,
                               resource_pool=resource_pool_id,
                               storage_path_prefix=storage_path_prefix)
    task = submit_request(
        api_instance=api_instance,
        request=task_request,
    )

    task_id = task["id"]
    log_task_info(task_id, method_name, params, resource_pool)

    with blocking_task_context(api_instance, task_id):
        if task["status"] == "pending-input":

            upload_input(
                api_instance=api_instance,
                original_params=params,
                task_id=task_id,
                type_annotations=type_annotations,
            )

    return task_id


def invoke_async_api(method_name: str,
                     params,
                     type_annotations: Dict[Any, Type],
                     resource_pool: Optional[
                         types.ComputationalResources] = None,
                     storage_path_prefix: Optional[str] = "") -> str:
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

    Return:
        Returns the task id.
    """

    api_config = validate_api_key(inductiva.api_key)

    request_params = get_validate_request_params(
        original_params=params,
        type_annotations=type_annotations,
    )

    with ApiClient(api_config) as client:
        api_instance = TasksApi(client)

        task_id = submit_task(api_instance=api_instance,
                              method_name=method_name,
                              request_params=request_params,
                              resource_pool=resource_pool,
                              storage_path_prefix=storage_path_prefix,
                              params=params,
                              type_annotations=type_annotations)

    return task_id


def compare_client_and_backend_versions(client_version: str):
    """ Compares the provided client version 7with the backend API version.

    Sends a GET request to the backend API's version comparison endpoint 
    with the client version as a parameter. Evaluates the response to 
    determine if the client version is compatible with the backend version. 
    Raises exceptions for communication issues or incompatibility.

    Parameters:
    - client_version (str): The version of the client to be compared with the 
                            backend version.

    Raises:
    - RuntimeError: If the API cannot be reached, or if the client version is 
      incompatible with the backend version, or for other general failures.
    """
    api_config = Configuration(host=inductiva.api_url)

    with ApiClient(api_config) as client:
        api_instance = VersionApi(client)
        query_params = {"client_version": client_version}

        try:
            api_instance.compare_client_and_backend_versions(
                query_params=query_params)

        except (MaxRetryError, NewConnectionError) as exc:
            raise RuntimeError(
                "Failed to reach the API. "
                "Please check your connection and try again.") from exc

        except ApiException as e:
            if e.status == 406:
                raise RuntimeError(
                    f"Client version {client_version} is not compatible "
                    f"with API version {e.headers['version']}.\n"
                    "Please update the client version.") from e
            raise RuntimeError(e) from e

        except Exception as e:
            raise RuntimeError(
                f"Failed to compare client and API versions. {e}") from e
