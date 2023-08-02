"""Functions for running simulations via Inductiva Web API."""
import pathlib
from typing import Any, Optional, Union
from uuid import UUID

from inductiva.tasks import Task
from inductiva.api.methods import invoke_async_api


def run_simulation(
    api_method_name: str,
    input_dir: pathlib.Path,
    output_dir: Optional[pathlib.Path] = None,
    resource_pool_id: Optional[UUID] = None,
    run_async: bool = False,
    **kwargs: Any,
) -> Union[pathlib.Path, Task]:
    """Run a simulation via Inductiva Web API."""

    params = {
        "sim_dir": input_dir,
        **kwargs,
    }
    type_annotations = {
        "sim_dir": pathlib.Path,
    }

    task_id = invoke_async_api(api_method_name,
                               params,
                               type_annotations,
                               resource_pool_id=resource_pool_id)
    task = Task(task_id)
    if not isinstance(task_id, str):
        raise RuntimeError(
            f"Expected result to be a string with task_id, got {type(task_id)}")
    
    # Blocking call for sync execution
    if not(run_async):
        with task:
            task.wait()

    return task
