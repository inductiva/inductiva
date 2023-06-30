"""Methods that interact with the lower-level inductiva-web-api-client.

The relevant function for usage outside of this file is invoke_api().
Check the demos directory for examples on how it is used.
"""
import os
import pathlib
import signal
import time
from contextlib import contextmanager
from typing import Any, Dict, Optional, Type
from uuid import UUID

from absl import logging

import inductiva
from inductiva.client import ApiClient, ApiException, Configuration
from inductiva.client.apis.tags.tasks_api import TasksApi
from inductiva.client.models import (BodyUploadTaskInput, TaskRequest,
                                     TaskStatus)
from inductiva.exceptions import RemoteExecutionError
from inductiva.types import Path
from inductiva.utils.data import (extract_output, get_validate_request_params,
                                  pack_input, unpack_output)
from inductiva.utils.meta import get_method_name, get_type_annotations


def validate_api_key(api_key: Optional[str]) -> Configuration:
    """Validates the API key and returns API configuration"""
    if inductiva.api_key is None:
        raise ValueError(
            "No API Key specified. "
            "Set it in the code with \"inductiva.api_key = <YOUR_SECRET_KEY>\""
            " or set the INDUCTIVA_API_KEY environment variable.")

    api_config = Configuration(host=inductiva.api_url)
    api_config.api_key["APIKeyHeader"] = api_key

    return api_config


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
    input_zip_path = pack_input(
        params=original_params,
        type_annotations=type_annotations,
        zip_name=task_id,
    )

    logging.info("Uploading input files...")
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

    os.remove(input_zip_path)


def download_output(api_instance: TasksApi,
                    task_id,
                    output_dir: Optional[Path] = None):
    """Downloads the output of a given task from the API.

    Args:
        api_instance: Instance of TasksApi used to send necessary requests.
        task_id: ID of the task.

    Return:
        Downloads and extracts the data to the client.
    """
    logging.info("Downloading output...")
    try:
        api_response = api_instance.download_task_output(
            path_params={"task_id": task_id},
            stream=True,
        )
        logging.info("Output downloaded.")
    except ApiException as e:
        raise e

    logging.debug("Downloaded output to %s", api_response.body.name)

    if output_dir is None:
        output_dir = os.path.join(inductiva.output_dir, task_id)

    logging.info("Extracting output ZIP file to \"%s\"...", output_dir)
    result_list = extract_output(api_response.body.name, output_dir)
    logging.info("Output extracted.")

    return result_list


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


def invoke_api_from_fn_ptr(
    params,
    function_ptr,
    output_dir: Optional[Path] = None,
    resource_pool_id: Optional[UUID] = None,
):
    """Perform a task remotely defined by a function pointer."""
    type_annotations = get_type_annotations(function_ptr)
    method_name = get_method_name(function_ptr)
    return invoke_api(
        method_name,
        params,
        output_dir=output_dir,
        type_annotations=type_annotations,
        return_type=type_annotations["return"],
        resource_pool_id=resource_pool_id,
    )


def submit_task(api_instance, task_request, params, type_annotations):
    """Submit a task and send input files to the API."""

    task = submit_request(
        api_instance=api_instance,
        request=task_request,
    )

    task_id = task["id"]
    logging.info("This simulation has the task_id: %s", task_id)

    with blocking_task_context(api_instance, task_id):
        if task["status"] == "pending-input":

            upload_input(
                api_instance=api_instance,
                original_params=params,
                task_id=task_id,
                type_annotations=type_annotations,
            )

            logging.info("Request submitted.")

    return task_id


def invoke_api(method_name: str,
               params,
               type_annotations: Dict[Any, Type],
               output_dir: Optional[Path] = None,
               return_type: Type = pathlib.Path,
               resource_pool_id: Optional[UUID] = None):
    """Perform a task remotely via Inductiva's Web API.

    Currently, the implementation handles the whole flow of the task execution,
    and blocks until the task finishes execution.
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
        4. Block execution, polling using the "GET task/{task_id}/status"
            endpoint until the returned status is "success".
        5. Download output ZIP with "GET task/{task_id}/output".
        6. Unpack output ZIP and return correct value based on the
            type annotations of `function_ptr`.

    Args:
        request: Request sent to the API for validation.
        input_dir: Directory containing the input files to be uploaded.
        output_dir: Directory where to place the output files.
        return_type: Type of the return value of the task, for unpacking.

    Return:
        Returns the output of the task.
    """
    api_config = validate_api_key(inductiva.api_key)

    request_params = get_validate_request_params(
        original_params=params,
        type_annotations=type_annotations,
    )

    task_request = TaskRequest(
        method=method_name,
        params=request_params,
        resource_pool=resource_pool_id,
    )

    with ApiClient(api_config) as client:
        api_instance = TasksApi(client)

        task_id = submit_task(api_instance, task_request, params,
                              type_annotations)

        block_until_status_is(
            api_instance=api_instance,
            task_id=task_id,
            desired_status={"started"},
        )
        logging.info("An executer has picked up the request")
        logging.info("The requested task is being executed remotely")

        # While the task is executing, use a context manager that kills the
        # task if some exception or SIGINT is caught.
        with blocking_task_context(api_instance, task_id):
            status = block_until_finish(
                api_instance=api_instance,
                task_id=task_id,
            )

            if status == "success":
                logging.info("Task executed successfuly.")
            else:
                logging.info("Task failed.")

        result_list = download_output(api_instance, task_id, output_dir)

        if status == "failed":
            raise RemoteExecutionError(f"""Remote execution failed.
        Find more details in: \"{os.path.abspath(output_dir)}\".""")

    return unpack_output(result_list, output_dir, return_type)


def invoke_async_api(method_name: str,
                     params,
                     type_annotations: Dict[Any, Type],
                     resource_pool_id: Optional[UUID] = None):
    """Perform a task asyc and remotely via Inductiva's Web API.

    Submits a simulation async to the API and returns the task id.
    The flow is similar to invoke_api, but it does not block until the task
    is finished.
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

    task_request = TaskRequest(
        method=method_name,
        params=request_params,
        resource_pool=resource_pool_id,
    )

    with ApiClient(api_config) as client:
        api_instance = TasksApi(client)

        task_id = submit_task(api_instance=api_instance,
                              task_request=task_request,
                              params=params,
                              type_annotations=type_annotations)

    return task_id


def run_simulation(
    api_method_name: str,
    input_dir: pathlib.Path,
    output_dir: pathlib.Path,
    params: Dict[str, Any],
    resource_pool_id: Optional[UUID] = None,
) -> pathlib.Path:
    """Run a simulation synchronously via Inductiva Web API."""

    params = {
        "sim_dir": input_dir,
        **params,
    }
    type_annotations = {
        "sim_dir": pathlib.Path,
    }

    result = invoke_api(
        api_method_name,
        params,
        type_annotations,
        output_dir=output_dir,
        return_type=pathlib.Path,
        resource_pool_id=resource_pool_id,
    )

    if not isinstance(result, pathlib.Path):
        raise RuntimeError(f"Expected result to be a Path, got {type(result)}")

    return result


def run_async_simulation(
    api_method_name: str,
    input_dir: pathlib.Path,
    params: Dict[str, Any],
    resource_pool_id: Optional[UUID] = None,
) -> str:
    """Run a simulation asynchronously via Inductiva Web API."""

    params = {
        "sim_dir": input_dir,
        **params,
    }
    type_annotations = {
        "sim_dir": pathlib.Path,
    }

    result = invoke_async_api(api_method_name,
                              params,
                              type_annotations,
                              resource_pool_id=resource_pool_id)

    if not isinstance(result, str):
        raise RuntimeError(
            f"Expected result to be a string with task_id, got {type(result)}")

    return result
