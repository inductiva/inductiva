"""Functions for running simulations via Inductiva Web API."""
import pathlib
from typing import Any, List, Optional

import logging

from inductiva import tasks, types, utils
from inductiva.api import methods


def run_simulation(
    simulator: str,
    input_dir: Optional[pathlib.Path],
    *,
    machine_group: types.ComputationalResources,
    resubmit_on_preemption: bool = False,
    storage_dir: Optional[str] = "",
    simulator_obj=None,
    remote_assets: Optional[List[str]] = None,
    simulator_name_alias: Optional[str] = None,
    project_name: Optional[str] = None,
    **kwargs: Any,
) -> tasks.Task:
    """Run a simulation via Inductiva Web API."""

    params = {
        "sim_dir": utils.data.INPUT_DIRNAME,
        **kwargs,
    }

    if not remote_assets:
        remote_assets = []

    container_image = kwargs.get("container_image", None)

    if (machine_group.allow_auto_start and not machine_group.started):
        logging.info("\n■ The computational resource is not started."
                     " Starting it now.\n")
        machine_group.start()

    task_id = methods.submit_task(simulator,
                                  input_dir=input_dir,
                                  machine_group=machine_group,
                                  params=params,
                                  storage_path_prefix=storage_dir,
                                  resubmit_on_preemption=resubmit_on_preemption,
                                  container_image=container_image,
                                  simulator_name_alias=simulator_name_alias,
                                  simulator_obj=simulator_obj,
                                  remote_assets=remote_assets,
                                  project_name=project_name)

    logging.info("■ Task %s submitted to the queue of the %s.", task_id,
                 machine_group)

    if not isinstance(task_id, str):
        raise RuntimeError(
            f"Expected result to be a string with task_id, got {type(task_id)}")
    task = tasks.Task(task_id)

    position = task.get_position_in_queue()
    if position is not None:
        pos_info = f"Number of tasks ahead in the queue: {position}."
    else:
        pos_info = f"Task {task_id} does not have queue information."

    logging.info(
        "%s\n"
        "· Consider tracking the status of the task via CLI:"
        "\n\tinductiva tasks list --id %s\n"
        "· Or, tracking the logs of the task via CLI:"
        "\n\tinductiva logs %s\n"
        "· Or, track the task files in real time with:"
        "\n\tinductiva tasks list-files %s\n"
        "· You can also get more information "
        "about the task via the CLI command:"
        "\n\tinductiva tasks info %s\n\n", pos_info, task_id, task_id, task_id,
        task_id)

    return task
