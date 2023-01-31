"""
Methods that interact with the lower-level inductiva-web-api-client.
"""
import os
import time
from absl import logging

from inductiva_web_api_client import ApiClient, ApiException
from inductiva_web_api_client.apis.tags.tasks_api import TasksApi
from inductiva_web_api_client.models import TaskRequest

from inductiva.utils.data import get_validate_request_params, pack_input, unpack_output
from inductiva.utils.meta import get_type_annotations, get_method_name
from inductiva.config import Configuration

configuration = None


def init(address, output_dir):
    global configuration
    configuration = Configuration(address=address, output_dir=output_dir)


def is_initialized():
    return configuration is not None


def block_until_finish(api_instance, task_id, sleep_secs=0.5):
    while True:
        try:
            api_response = \
                api_instance.get_task_status_task_task_id_status_get(
                    path_params={"task_id": task_id},
                    )
            logging.debug(api_response)
        except ApiException as e:
            raise e

        # If status is success, then stop polling
        if api_response.body["status"] == "success":
            break

        time.sleep(sleep_secs)

    return api_response


def invoke_api(params, function_ptr):
    if not is_initialized():
        raise Exception("Inductiva not initialized.")

    type_annotations = get_type_annotations(function_ptr)

    request_params = get_validate_request_params(
        original_params=params,
        type_annotations=type_annotations,
    )

    with ApiClient(configuration.api_config) as client:
        api_instance = TasksApi(client)

        task_request = TaskRequest(
            method=get_method_name(function_ptr),
            params=request_params,
        )

        # Submit task
        try:
            api_response = api_instance.submit_task_task_submit_post(
                body=task_request,)
            logging.debug(api_response)
        except ApiException as e:
            logging.exception(
                "Exception when calling TasksApi->submit_task: %s", e)
            raise e

        task_id = api_response.body["id"]

        # Upload input as zip if required
        if api_response.body["status"] == "pending-input":
            input_zip_path = pack_input(
                params=params,
                type_annotations=type_annotations,
                zip_name=task_id,
            )

            logging.debug("Uploading input zip ...")
            try:
                with open(input_zip_path, "rb") as zip_fp:
                    api_response = \
                        api_instance.upload_task_input_task_task_id_input_post(
                            path_params=dict(task_id=task_id),
                            body=dict(file=zip_fp),
                        )
                    logging.info(api_response)
            except ApiException as e:
                logging.exception(
                    "Exception when calling TasksApi->upload_task_inputs: %s",
                    e)
                raise e

            os.remove(input_zip_path)

        response = block_until_finish(
            api_instance=api_instance,
            task_id=task_id,
        )

        logging.debug("Request status: %s", response.body["status"])

        # Download output
        try:
            api_response = api_instance.get_task_output_task_task_id_output_get(
                path_params=dict(task_id=task_id),
                stream=True,
            )
            logging.debug(api_response)
        except ApiException as e:
            raise e

    logging.debug("Downloaded output to %s", api_response.body.name)

    return unpack_output(
        zip_path=api_response.body.name,
        output_path=os.path.join(configuration.output_dir, task_id),
        return_type=type_annotations["return"],
    )
