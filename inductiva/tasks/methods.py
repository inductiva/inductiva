"""Methods to interact with the tasks submitted to the API."""
import pathlib
from typing import Optional, Type

from absl import logging

import inductiva
from inductiva import api
from inductiva.types import Path
from inductiva.client import ApiClient
from inductiva.client.apis.tags.tasks_api import TasksApi
from inductiva.utils.data import (unpack_output)


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
