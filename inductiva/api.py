"""Methods that interact with the lower-level inductiva-web-api-client.

The most relevant functions that for usage outside of this file
are init() and invoke_api(). Check the demos directory and the `math.py`
file for examples on how they are used.
"""
import os
import time
from absl import logging

from inductiva_web_api_client import ApiClient, ApiException
from inductiva_web_api_client.apis.tags.tasks_api import TasksApi
from inductiva_web_api_client.models import TaskRequest, TaskStatus

from inductiva.utils.data import get_validate_request_params, pack_input, unpack_output
from inductiva.utils.meta import get_type_annotations, get_method_name
from inductiva.config import Configuration

DEFAULT_OUTPUT_DIR = "inductiva_output"
DEFAULT_API_URL = "http://api.inductiva.ai"

configuration = Configuration(address=DEFAULT_API_URL,
                              output_dir=DEFAULT_OUTPUT_DIR)


def update_config(address=DEFAULT_API_URL, output_dir=DEFAULT_OUTPUT_DIR):
    """Update the configuration.

    Args:
        address: Address (including port) where to connect to the API.
        output_dir: Path in which to store outputs of the executed tasks.
            Outputs of a given task will be stored in a child directory of
            `output_dir` with the name equal to the ID of the task.
    """
    global configuration
    configuration = Configuration(address=address, output_dir=output_dir)


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


def block_until_finish(api_instance,
                       task_id: str,
                       sleep_secs=0.5) -> TaskRequest:
    """Block until a task executing remotely finishes execution.

    Args:
        api_instance: Instance of TasksApi used to send necessary requests.
        task_id: ID of the task to wait for.
        sleep_secs: Time in secs between polling requests. Defaults to 0.5.

    Return
        Returns info related to the task, containing two fields,
        "id" and "status".
    """

    logging.debug("Blocking until task is finished ...")
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
        if api_response.body["status"] == "success":
            break

        time.sleep(sleep_secs)

    return api_response.body


def invoke_api(params, function_ptr):
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

    with ApiClient(configuration.api_config) as client:
        api_instance = TasksApi(client)

        task = submit_request(
            api_instance=api_instance,
            original_params=params,
            function_ptr=function_ptr,
        )

        task_id = task["id"]

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

    return unpack_output(
        zip_path=output_zip_path,
        output_dir=os.path.join(configuration.output_dir, task_id),
        return_type=type_annotations["return"],
    )
