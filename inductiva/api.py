"""
Methods that interact with the lower-level inductiva-web-api-client.
"""
import os
import inspect
import numpy as np
import time
import json
import sys
import tempfile
import zipfile
import shutil
import time
from absl import logging
from pprint import pprint

from inductiva_web_api_client import Configuration
from inductiva_web_api_client import ApiClient, ApiException
from inductiva_web_api_client.apis.tags.tasks_api import TasksApi
from inductiva_web_api_client.models import TaskRequest
from inductiva_web_api_client.schemas import Unset

configuration = Configuration()


def init(address):
    global configuration
    configuration = Configuration(host=address)
    # TODO: `init` should probably perform additional validations,
    # such as some kind of "ping" to check if the address is valid.


def invoke_api(params, function_ptr):
    args_spec = inspect.getfullargspec(function_ptr)
    type_annotations = args_spec.annotations

    original_params = params.copy()
    for variable in params:
        param_type = type_annotations[variable]
        if param_type == np.ndarray:
            params[variable] = {
                "shape": params[variable].shape,
            }

    # NOTE: Doing everything inside the same context manager.
    # Is this the right approach?
    with ApiClient(configuration) as client:
        api_instance = TasksApi(client)

        # remove name of the main module (inductiva)
        # e.g. inductiva.math.sum becomes math.sum
        module_name = ".".join(function_ptr.__module__.split(".")[1:])

        task_request = TaskRequest(
            method=f"{module_name}.{function_ptr.__name__}",
            params=params,
        )

        try:
            # Submit task
            api_response = api_instance.submit_task_task_submit_post(
                body=task_request,)
            logging.debug(api_response)

        # TODO: Currently, the API uses error code 400 if the requested
        # method doesn't exist,
        # and error 403 if the client has no permission/credits
        # to execute the requested task
        # Those errors are caught as exceptions here
        # the exception is being thrown again for now
        except ApiException as e:
            logging.exception(
                "Exception when calling TasksApi->submit_task: %s", e)
            raise e

        # TODO: Response should have information if any large data
        # should be uploaded as file.
        task_id = api_response.body["id"]

        if api_response.body["status"] == "pending-input":
            with tempfile.TemporaryDirectory() as tmpdir_path:
                for variable in params:
                    param_type = type_annotations[variable]
                    if param_type == np.ndarray:
                        param_filename = f"{variable}.npy"
                        param_fullpath = os.path.join(tmpdir_path, param_filename)
                        np.save(param_fullpath, original_params[variable])
                        params[variable] = param_filename # specify path that contains the given argument

                        logging.info("Stored %s to %s", variable, param_fullpath)


                with open(os.path.join(tmpdir_path, "input.json"), "w", encoding="UTF-8") as fp:
                    json.dump(params, fp)

                zip_path = shutil.make_archive(os.path.join(tempfile.gettempdir(), f"{task_id}"), "zip", tmpdir_path)
                logging.info("Compressed inputs to %s", zip_path)

            logging.info("Uploading %s", zip_path)
            try:
                api_response = api_instance.upload_task_input_task_task_id_input_post(
                    path_params=dict(task_id=task_id),
                    body=dict(file=open(zip_path, "rb")),
                )
                logging.info(api_response)
            except ApiException as e:
                logging.exception(
                    "Exception when calling TasksApi->upload_task_inputs: %s", e)
                raise e

            os.remove(zip_path)

        # Poll status until it is "successful"
        # Is there a way to do this without polling?
        while True:
            try:
                api_response = \
                    api_instance.get_task_status_task_task_id_status_get(
                        path_params={"task_id": task_id},)
                logging.debug(api_response)
            except ApiException as e:
                raise e

            # If status is success, then stop polling
            if api_response.body["status"] == "success":
                break

            time.sleep(0.5)

        # Now that the task has been marked as "success", get output
        try:
            api_response = api_instance.get_task_output_task_task_id_output_get(
                path_params=dict(task_id=task_id), stream=True,
            )
            logging.debug(api_response)
        except ApiException as e:
            raise e

    with tempfile.TemporaryDirectory() as tmpdir_path:
        with open(os.path.join(tmpdir_path, 'output.zip'), 'wb') as fp:
            fp.write(api_response.response.data)

        with zipfile.ZipFile(os.path.join(tmpdir_path, 'output.zip'), "r") as zip_fp:
            zip_fp.extractall(os.path.join(tmpdir_path, 'output'))

        with open(os.path.join(tmpdir_path, 'output', 'output.json'), "r") as fp:
            result_list = json.load(fp)

        logging.info("Output dict %s", result_list)

        # TODO: tratar de multiplos returns
        return_type = type_annotations["return"]
        if return_type == np.ndarray:
            return_value = np.load(os.path.join(tmpdir_path, 'output', result_list[0]))
        else:
            return_value = return_type(result_list[0])

        return return_value

