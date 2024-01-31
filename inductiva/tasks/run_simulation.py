"""Functions for running simulations via Inductiva Web API."""
import os

import pathlib
from typing import Any, Optional
import json
import threading

from absl import logging

from inductiva import tasks, types
from inductiva.api import methods
from inductiva.utils import format_utils, files

TASK_METADATA_FILENAME = "task_metadata.json"

_metadata_lock = threading.RLock()


def run_simulation(
    api_method_name: str,
    input_dir: pathlib.Path,
    computational_resources: Optional[types.ComputationalResources] = None,
    storage_dir: Optional[types.Path] = "",
    api_invoker=None,
    extra_metadata=None,
    **kwargs: Any,
) -> tasks.Task:
    """Run a simulation via Inductiva Web API."""

    params = {
        "sim_dir": input_dir,
        **kwargs,
    }
    type_annotations = {
        "sim_dir": pathlib.Path,
    }

    if api_invoker is None:
        api_invoker = methods.invoke_async_api

    task_id = api_invoker(
        api_method_name,
        params,
        type_annotations,
        resource_pool=computational_resources,
        storage_path_prefix=storage_dir,
    )

    if computational_resources is not None:
        logging.info("Task %s submitted to the queue of the %s.", task_id,
                     computational_resources)
    else:
        logging.info("Task %s submitted to the default queue.", task_id)

    task = tasks.Task(task_id)
    if not isinstance(task_id, str):
        raise RuntimeError(
            f"Expected result to be a string with task_id, got {type(task_id)}")

    if not format_utils.getenv_bool("DISABLE_TASK_METADATA_LOGGING", False):
        machine_group_id = None
        if computational_resources is not None:
            machine_group_id = computational_resources.id

        metadata = {
            "api_method_name": api_method_name.split(".")[1],
            "machine_group_id": machine_group_id,
            "storage_dir": str(storage_dir),
            **kwargs,
        }
        if extra_metadata is not None:
            metadata = {**metadata, **extra_metadata}

        with _metadata_lock:
            _save_metadata({
                **{
                    "task_id": task_id,
                    "input_dir": str(input_dir)
                },
                **metadata
            })
        logging.info(
            "Task %s configurations metadata saved to the tasks metadata file "
            "%s in the current working directory.", task_id,
            TASK_METADATA_FILENAME)

    logging.info(
        "Consider tracking the status of the task via CLI:"
        "\n\tinductiva tasks list --task-id %s", task_id)
    logging.info(
        "Or, tracking the logs of the task via CLI:"
        "\n\tinductiva logs %s", task_id)

    return task


def _save_metadata(metadata, mode="a"):
    """Appends metadata to the TASK_METADATA_FILENAME in the cwd."""

    file_path = files.resolve_output_path(TASK_METADATA_FILENAME)
    if not os.path.exists(file_path.parent):
        os.mkdir(file_path.parent)
    with open(file_path, mode, encoding="utf-8") as f:
        json.dump(metadata, f)
        f.write("\n")
    logging.info("Simulation metadata logged to: %s", file_path)
