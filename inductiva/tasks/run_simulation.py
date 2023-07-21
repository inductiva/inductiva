"""Functions for running simulations via Inductiva Web API."""
import pathlib
from typing import Any, Dict, Optional
from uuid import UUID

from inductiva.tasks import Task
from inductiva.api.methods import invoke_async_api


def run_simulation(
    api_method_name: str,
    input_dir: pathlib.Path,
    output_dir: pathlib.Path,
    params: Dict[str, Any],
    resource_pool_id: Optional[UUID] = None,
) -> pathlib.Path:
    """Run a simulation synchronously via Inductiva Web API."""
    task = run_async_simulation(api_method_name, input_dir, params,
                                resource_pool_id)
    with task:
        task.wait()

    return task.download_output(output_dir)


def run_async_simulation(
    api_method_name: str,
    input_dir: pathlib.Path,
    params: Dict[str, Any],
    resource_pool_id: Optional[UUID] = None,
) -> Task:
    """Run a simulation asynchronously via Inductiva Web API."""

    params = {
        "sim_dir": input_dir,
        **params,
    }
    type_annotations = {
        "sim_dir": pathlib.Path,
    }

    task_id = invoke_async_api(api_method_name,
                               params,
                               type_annotations,
                               resource_pool_id=resource_pool_id)

    if not isinstance(task_id, str):
        raise RuntimeError(
            f"Expected result to be a string with task_id, got {type(task_id)}")

    return Task(task_id)
