"""Base class for low-level simulators."""
import os
import pathlib
from typing import Optional, Type

from absl import logging

import inductiva
from inductiva import api
from inductiva.types import Path
from inductiva.client import ApiClient, Configuration
from inductiva.client.apis.tags.tasks_api import TasksApi
from inductiva.utils.data import (extract_output, unpack_output)


def check_status(task_id: str):
    """Check the status of a simulation."""

    if inductiva.api_key is None:
        raise ValueError(
            "No API Key specified. "
            "Set it in the code with \"inductiva.api_key = <YOUR_SECRET_KEY>\""
            " or set the INDUCTIVA_API_KEY environment variable.")

    api_config = Configuration(host=inductiva.api_url)
    api_config.api_key["APIKeyHeader"] = inductiva.api_key
    logging.info("Checking status of task %s", task_id)

    with ApiClient(api_config) as client:
        api_instance = TasksApi(client)

        status = api.methods.check_status(api_instance=api_instance,
                                          task_id=task_id)

    return status


def fetch_task(task_id: str,
               output_dir: Optional[Path] = None,
               return_type: Type = pathlib.Path) -> pathlib.Path:
    """Fetch the results of a simulation."""

    if inductiva.api_key is None:
        raise ValueError(
            "No API Key specified. "
            "Set it in the code with \"inductiva.api_key = <YOUR_SECRET_KEY>\""
            " or set the INDUCTIVA_API_KEY environment variable.")

    api_config = Configuration(host=inductiva.api_url)
    api_config.api_key["APIKeyHeader"] = inductiva.api_key
    logging.info("Connecting to API")

    with ApiClient(api_config) as client:
        api_instance = TasksApi(client)

        status = api.check_status(api_instance=api_instance, task_id=task_id)

        if status == "started":
            logging.info("The task is still being executed.")
        elif status == "submitted":
            logging.info("The task is waiting for resources...")
        else:

            logging.info("Downloading output...")
            output_zip_path = api.download_output(
                api_instance=api_instance,
                task_id=task_id,
            )
            logging.info("Output downloaded.")

            if output_dir is None:
                output_dir = os.path.join(inductiva.output_dir, task_id)

            logging.info("Extracting output ZIP file to \"%s\"...", output_dir)
            result_list = extract_output(output_zip_path, output_dir)
            logging.info("Output extracted.")

            return unpack_output(result_list, output_dir, return_type)
