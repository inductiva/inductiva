import inspect
import numpy as np
import time
from absl import logging

from inductiva_web_api_client import Configuration
from inductiva_web_api_client import ApiClient, ApiException
from inductiva_web_api_client.apis.tags.tasks_api import TasksApi
from inductiva_web_api_client.models import TaskRequest

configuration = Configuration()


def init(address):
    global configuration
    configuration = Configuration(host=address)
    # Note: `init` should probably perform additional validations,
    # such as some kind of "ping" to check if the address is valid.


def invoke_api(params, function_ptr):
    args_spec = inspect.getfullargspec(function_ptr)
    type_annotations = args_spec.annotations

    original_params = params.copy()
    for variable in params:
        param_type = type_annotations[variable]
        if param_type == np.ndarray:
            params[variable] = {"shape": params[variable].shape}

    # NOTE: Doing everything inside the same context manager.
    # Is this the right approach?
    with ApiClient(configuration) as client:
        api_instance = TasksApi(client)

        task_request = TaskRequest(
            method=".".join(function_ptr.__module__.split(".")[1:]) + "." +
            function_ptr.__name__,  # removes inductiva. from the module name
            params=params,
        )

        try:
            # Submit task
            api_response = api_instance.submit_task_task_submit_post(
                body=task_request,)
            logging.debug(api_response)

        # Currently, the API uses error code 400 if the method requested doesn't exist,
        # and error 403 if the client has no permission/credits to execute the requested task
        # Those errors are caught as execptions here
        # the exception is being thrown again for now
        except ApiException as e:
            logging.exception(
                "Exception when calling TasksApi->submit_task_task_submit_post: %s",
                e)
            raise e

        # Response should have information if any large data should be uploaded as file.
        # That's WIP; right it is assumed that the request contains everything needed to perform the task
        task_id = api_response.body["id"]

        # Poll status until it is "successful"
        while True:
            try:
                api_response = api_instance.get_task_status_task_task_id_status_get(
                    path_params={"task_id": task_id},)
                logging.debug(api_response)
            except ApiException as e:
                raise e

            # If status is success, then stop polling
            if api_response.body["status"] == "success":
                break

            time.sleep(0.5)

        try:
            api_response = api_instance.get_task_output_task_task_id_output_get(
                path_params={"task_id": task_id},)
            logging.debug(api_response)
        except ApiException as e:
            raise e

    # So far experimented only with the dummy.sum endpoint: it specifies the result in a `c` field.
    # All endpoints that return a single result should use the same field (e.g. `result`).
    # Need to figure out how to handle endpoints that generate multiple outputs.
    result = api_response.body["c"]

    return_type = type_annotations["return"]
    if return_type == np.ndarray:
        return_value = np.frombuffer(result)
    else:
        return_value = return_type(result)
    return return_value
