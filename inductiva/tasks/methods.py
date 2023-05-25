"""Methods to interact with the tasks submitted to the API."""
import pathlib
from typing import Dict, List, Optional, Type

from absl import logging

import inductiva
from inductiva import api
from inductiva.client import ApiClient, ApiException
from inductiva.client.apis.tags.tasks_api import TasksApi
from inductiva.client.model.task_status_code import TaskStatusCode
from inductiva.types import Path
from inductiva.utils.data import unpack_output


def get_task_status(task_id: str):
    """Get the status of a simulation."""

    api_config = api.validate_api_key(inductiva.api_key)

    with ApiClient(api_config) as client:
        api_instance = TasksApi(client)

        status = api.get_task_status(api_instance, task_id)
    logging.info("Task status: %s", status)

    return status


def fetch_task_output(task_id: str,
                      output_dir: Optional[Path] = None,
                      return_type: Type = pathlib.Path) -> Type:
    """Fetch the results of a simulation."""
    logging.info("Fetching the simulation output for task %s", task_id)
    api_config = api.validate_api_key(inductiva.api_key)

    with ApiClient(api_config) as client:
        api_instance = TasksApi(client)

        status = api.get_task_status(api_instance, task_id)

        logging.info(status)
        if status == "started":
            logging.info("The task is still being executed.")
            return_type = None
        elif status == "submitted":
            logging.info("The task is waiting for resources...")
            return_type = None
        else:
            result_list = api.download_output(api_instance, task_id, output_dir)

            return unpack_output(result_list, output_dir, return_type)


def get_task_info(task_id) -> Dict:
    """Get info about a task by its ID."""
    api_config = api.validate_api_key(inductiva.api_key)

    with ApiClient(api_config) as client:
        api_instance = TasksApi(client)

        try:
            # Get User Tasks
            api_response = api_instance.get_task(
                path_params={"task_id": task_id})

            # Convert DynamicSchema to dict
            return {**api_response.body}

        except ApiException as e:
            raise e


def get_tasks_info(status: Optional[str] = None,
                   page=1,
                   per_page=10) -> List[Dict]:
    """Get information about a user's tasks on the API.

    Tags can be filtered by a status. Results are paginated indexed from 1.
    """
    api_config = api.validate_api_key(inductiva.api_key)

    with ApiClient(api_config) as client:
        api_instance = TasksApi(client)

        query_params = {
            "page": page,
            "per_page": per_page,
        }

        if status is not None:
            query_params["status"] = TaskStatusCode(status)

        try:
            # Get User Tasks
            api_response = api_instance.get_user_tasks(
                query_params=query_params,)

            return [{**task} for task in api_response.body]

        except ApiException as e:
            raise e
