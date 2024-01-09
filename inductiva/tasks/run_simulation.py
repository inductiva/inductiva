"""Functions for running simulations via Inductiva Web API."""
import pathlib
from typing import Any, Optional
import json

from absl import logging

from inductiva import tasks, resources, types
from inductiva.api import methods
from inductiva.utils import format_utils


def run_simulation(
    api_method_name: str,
    input_dir: pathlib.Path,
    machine_group: Optional[resources.MachineGroup] = None,
    storage_dir: Optional[types.Path] = "",
    save_sim_metadata_dir: Optional[types.Path] = None,
    extra_sim_metadata_to_save: Optional[dict] = None,
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

    if format_utils.getenv_bool("SAVE_SIM_METADATA", True):
        metadata = {
            **{
                "machine_group_id": resource_pool_id,
                "task_id": task_id,
                "input_dir": str(input_dir),
                "storage_dir": str(storage_dir),
            },
            **kwargs,
        }
        if extra_sim_metadata_to_save is not None:
            metadata = {**metadata, **extra_sim_metadata_to_save}
        _save_metadata(save_sim_metadata_dir, metadata)

    return task


def _save_metadata(save_sim_metadata_dir, metadata):
    """Save metadata to a .json file.

    Args:
      save_sim_metadata_dir: Directory where to save the metadata. If
      `None` it saves on the current working directory.
      metadata: Dictionary with the metadata to be saved.

    """
    if save_sim_metadata_dir is None:
        save_sim_metadata_dir = pathlib.Path().cwd()
    # Ensure that it is a pthlib object in case the user passed a string.
    save_sim_metadata_dir = pathlib.Path(save_sim_metadata_dir).resolve()

    task_id = metadata["task_id"]
    save_sim_metadata_dir /= f"{task_id}"
    save_sim_metadata_dir.mkdir(parents=True, exist_ok=True)
    file_path = save_sim_metadata_dir / "sim_metadata.json"

    with open(file_path, "w", encoding="utf-8") as f:
        json.dump(metadata, f)

    logging.info("Simulation metadata logged to: %s", file_path)
