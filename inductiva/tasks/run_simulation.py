"""Functions for running simulations via Inductiva Web API."""
import pathlib
from typing import Any, Optional
import json
import threading

from absl import logging

from inductiva import tasks, types
from inductiva.api import methods
from inductiva.utils import format_utils, files

TASK_METADATA_FILENAME = "task_metadata.json"
TASK_METADATA_FILENAME_UPLOAD = "uploaded_metadata.json"

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

    resource_pool_id = None
    if computational_resources is not None:
        resource_pool_id = computational_resources.id

    if api_invoker is None:
        api_invoker = methods.invoke_async_api

    if not format_utils.getenv_bool("DISABLE_TASK_METADATA_LOGGING", False):
        metadata = {
            "api_method_name": api_method_name.split(".")[1],
            "machine_group_id": resource_pool_id,
            "storage_dir": str(storage_dir),
            **kwargs,
        }
        if extra_metadata is not None:
            metadata = {**metadata, **extra_metadata}

        with _metadata_lock:
            _save_metadata(metadata,
                           mode="w",
                           path=pathlib.Path(input_dir) /
                           TASK_METADATA_FILENAME_UPLOAD)

    task_id = api_invoker(
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
        with _metadata_lock:
            _save_metadata({
                **{
                    "task_id": task_id,
                    "input_dir": str(input_dir)
                },
                **metadata
            })

    return task


def _save_metadata(metadata, mode="a", path=None):
    """Appends metadata to the TASK_METADATA_FILENAME in the cwd."""
    if path is None:
        file_path = files.resolve_path(TASK_METADATA_FILENAME)
    else:
        file_path = files.resolve_path(path)
    with open(file_path, mode, encoding="utf-8") as f:
        json.dump(metadata, f)
        f.write("\n")
    logging.info("Simulation metadata logged to: %s", file_path)
