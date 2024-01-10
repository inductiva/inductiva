"""Functions for running simulations via Inductiva Web API."""
import pathlib
from typing import Any, Optional
import json

from absl import logging

from inductiva import tasks, resources, types
from inductiva.api import methods
from inductiva.utils import format_utils

TASK_METADATA_FILENAME = "task_metadata.json"


def run_simulation(
    api_method_name: str,
    input_dir: pathlib.Path,
    machine_group: Optional[resources.MachineGroup] = None,
    storage_dir: Optional[types.Path] = "",
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

    resource_pool_id = None
    if machine_group is not None:
        resource_pool_id = machine_group.id

    task_id = methods.invoke_async_api(
        api_method_name,
        params,
        type_annotations,
        resource_pool_id=resource_pool_id,
        storage_path_prefix=storage_dir,
    )
    task = tasks.Task(task_id)
    if not isinstance(task_id, str):
        raise RuntimeError(
            f"Expected result to be a string with task_id, got {type(task_id)}")

    if not format_utils.getenv_bool("DISABLE_TASK_METADATA_LOGGING", False):
        metadata = {
            **{
                "task_id": task_id,
                "input_dir": str(input_dir),
                "storage_dir": str(storage_dir),
                "api_method_name": api_method_name,
                "machine_group_id": resource_pool_id,
            },
            **kwargs,
        }
        _save_metadata(metadata)

    return task


def _save_metadata(metadata):
    """Appends metadata to the TASK_METADATA_FILENAME in the cwd."""
    file_path = pathlib.Path().cwd() / TASK_METADATA_FILENAME
    with open(file_path, "a", encoding="utf-8") as f:
        json.dump(metadata, f)
        f.write("\n")
    logging.info("Simulation metadata logged to: %s", file_path)
