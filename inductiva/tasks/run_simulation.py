"""Functions for running simulations via Inductiva Web API."""
import pathlib
from typing import Any, Optional

from inductiva import tasks, types
from inductiva.api import methods


def run_simulation(
    api_method_name: str,
    input_dir: pathlib.Path,
    computational_resources: Optional[types.ComputationalResources] = None,
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
    if computational_resources is not None:
        resource_pool_id = computational_resources.id

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

    return task
