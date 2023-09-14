"""Functions for running simulations via Inductiva Web API."""
import pathlib
from typing import Any, Optional

from inductiva import tasks, resources
from inductiva.api import methods


def run_simulation(
    api_method_name: str,
    input_dir: pathlib.Path,
    machine_group: Optional[resources.MachineGroup] = None,
    run_async: bool = False,
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

    task_id = methods.invoke_async_api(api_method_name,
                                       params,
                                       type_annotations,
                                       resource_pool_id=resource_pool_id)
    task = tasks.Task(task_id)
    if not isinstance(task_id, str):
        raise RuntimeError(
            f"Expected result to be a string with task_id, got {type(task_id)}")

    # Blocking call for sync execution
    if not run_async:
        with task.sync_context():
            task.wait()

    return task
