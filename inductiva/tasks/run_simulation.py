"""Functions for running simulations via Inductiva Web API."""
import pathlib
from typing import Any, List, Optional

import logging

from inductiva import tasks, types
from inductiva.api import methods


def run_simulation(
    simulator: str,
    input_dir: Optional[pathlib.Path],
    *,
    computational_resources: types.ComputationalResources,
    resubmit_on_preemption: bool = False,
    storage_dir: Optional[str] = "",
    api_invoker=None,
    simulator_obj=None,
    input_resources: Optional[List[str]] = None,
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

    if not input_resources:
        input_resources = []

    container_image = kwargs.get("container_image", None)

    task_id = api_invoker(simulator,
                          params,
                          type_annotations,
                          computational_resources,
                          resubmit_on_preemption=resubmit_on_preemption,
                          container_image=container_image,
                          storage_path_prefix=storage_dir,
                          simulator_obj=simulator_obj,
                          input_resources=input_resources)
    logging.info("■ Task %s submitted to the queue of the %s.", task_id,
                 computational_resources)

    task = tasks.Task(task_id)
    if not isinstance(task_id, str):
        raise RuntimeError(
            f"Expected result to be a string with task_id, got {type(task_id)}")

    position = task.get_position_in_queue()
    if position is not None:
        position_message = f"Number of tasks ahead in the queue: {position}."
    else:
        position_message = f"Task {task_id} does not have queue information."

    logging.info(
        f"{position_message}\n"
        "· Consider tracking the status of the task via CLI:"
        f"\n\tinductiva tasks list --id {task_id}\n"
        "· Or, tracking the logs of the task via CLI:"
        f"\n\tinductiva logs {task_id}\n"
        "· You can also get more information "
        "about the task via the CLI command:"
        f"\n\tinductiva tasks info {task_id}\n\n")
    return task
