"""Methods that interact with the lower-level inductiva-web-api-client.

The relevant function for usage outside of this file is invoke_api().
Check the demos directory for examples on how it is used.
"""
import os
import signal
import time
from contextlib import contextmanager

from typing import Optional
from absl import logging
from inductiva_web_api_client import ApiClient, ApiException, Configuration
from inductiva_web_api_client.apis.tags.tasks_api import TasksApi
from inductiva_web_api_client.models import TaskRequest, TaskStatus
from inductiva.types import Path

import inductiva
from inductiva.utils.data import (get_validate_request_params, pack_input,
                                  unpack_output)
from inductiva.utils.meta import get_method_name, get_type_annotations


def submit_request(api_instance: TasksApi, original_params,
                   function_ptr) -> TaskStatus:
    """Submits a task request to the API.

    Args:
        api_instance: Instance of TasksApi used to send necessary requests.
        original_params: Params of the request passed by the user.
        function_ptr: Pointer to the function that defines the requested method.

    Return:
        Returns the body of the HTTP response.
        Contains two fields, "id" and "status".
    """
    request_params = get_validate_request_params(
        original_params=original_params,
        type_annotations=get_type_annotations(function_ptr),
    )

    task_request = TaskRequest(
        method=get_method_name(function_ptr),
        params=request_params,
    )

    try:
        api_response = api_instance.submit_task_task_submit_post(
            body=task_request)
    except ApiException as e:
        logging.exception("Exception when calling TasksApi->submit_task: %s", e)
        raise e

    logging.debug("Request status: %s", api_response.body["status"])

    return api_response.body


def upload_input(api_instance, task_id, original_params, type_annotations):
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

    logging.debug("Uploading input zip ...")
    try:
        with open(input_zip_path, "rb") as zip_fp:
            _ = api_instance.upload_task_input_task_task_id_input_post(
                path_params={"task_id": task_id},
                body={"file": zip_fp},
            )
    except ApiException as e:
        logging.exception(
            "Exception when calling TasksApi->upload_task_inputs: %s", e)
        raise e

    os.remove(input_zip_path)


def download_output(api_instance, task_id):
    """Downloads the output of a given task from the API.

    Args:
        api_instance: Instance of TasksApi used to send necessary requests.
        task_id: ID of the task.

    Return:
        Returns the path to the downloaded ZIP file.
    """
    try:
        api_response = api_instance.get_task_output_task_task_id_output_get(
            path_params={"task_id": task_id},
            stream=True,
        )
    except ApiException as e:
        raise e

    logging.debug("Downloaded output to %s", api_response.body.name)

    return api_response.body.name


def block_until_finish(api_instance, task_id: str) -> TaskRequest:
    """Block until a task executing remotely finishes execution.

    Args:
        api_instance: Instance of TasksApi used to send necessary requests.
        task_id: ID of the task to wait for.

    Returns:
        Returns info related to the task, containing two fields,
        "id" and "status".
    """
    logging.debug("Blocking until task is finished ...")
    return block_until_status_is(api_instance, task_id, "success")


def kill_task(api_instance, task_id: str):
    """Kill a task that is executing remotely.

    The function sends a Kill request to the API and then polls until the
    status of the task is killed.
    TODO: should polling be done here? Or should the API endpoint make sure
    that the task is killed instead of the client having to wait for
    confirmation?

    Args:
        api_instance: Instance of TasksApi used to send necessary requests.
        task_id: ID of the task to kill.

    Returns:
        Returns info related to the task, containing two fields,
        "id" and "status".
    """
    logging.debug("Sending kill task request ...")
    _ = api_instance.kill_task_task_task_id_kill_post(
        path_params={"task_id": task_id},)

    logging.debug("Blocking until task is killed ...")
    return block_until_status_is(api_instance, task_id, "killed")


def block_until_status_is(api_instance,
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
    while True:
        try:
            api_response = \
                api_instance.get_task_status_task_task_id_status_get(
                    path_params={"task_id": task_id},
                    )

            logging.debug("Request status: %s", api_response.body["status"])
        except ApiException as e:
            raise e

        # If status is success, then stop polling
        if api_response.body["status"] == desired_status:
            break

        time.sleep(sleep_secs)

    return api_response.body


@contextmanager
def blocking_task_context(api_instance, task_id):
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
        logging.debug("Caught exception: terminating blocking task ...")
        kill_task(api_instance, task_id)
        raise err
    except KeyboardInterrupt as err:
        logging.debug("Caught SIGINT: terminating blocking task ...")
        kill_task(api_instance, task_id)
        raise err
    finally:
        # Reset original SIGINT handler
        signal.signal(signal.SIGINT, original_sig)


def invoke_api(params, function_ptr, output_dir: Optional[Path] = None):
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
        params: Input params passed by the user.
        function_ptr: Pointer to the function that defines the method
            to be requested to the API.

    Return:
        Returns the output of the task.
    """

    type_annotations = get_type_annotations(function_ptr)

    api_config = Configuration(host=inductiva.api_url)

    with ApiClient(api_config) as client:
        api_instance = TasksApi(client)

        task = submit_request(
            api_instance=api_instance,
            original_params=params,
            function_ptr=function_ptr,
        )

        task_id = task["id"]

        # While the task is executing, use a context manager that kills the
        # task if some exception or SIGINT is caught.
        with blocking_task_context(api_instance, task_id):
            if task["status"] == "pending-input":
                upload_input(
                    api_instance=api_instance,
                    task_id=task_id,
                    original_params=params,
                    type_annotations=type_annotations,
                )

            _ = block_until_finish(
                api_instance=api_instance,
                task_id=task_id,
            )

        output_zip_path = download_output(
            api_instance=api_instance,
            task_id=task_id,
        )

    if output_dir is None:
        output_dir = os.path.join(inductiva.output_dir, task_id)

    return unpack_output(
        zip_path=output_zip_path,
        output_dir=output_dir,
        return_type=type_annotations["return"],
    )
