"""Methods to interact with the tasks submitted to the API."""
import json
from typing import Dict, List, Optional, Union, Sequence

import pandas as pd

import inductiva
from inductiva import api
from inductiva.client import ApiClient, ApiException
from inductiva.client.apis.tags.tasks_api import TasksApi
from inductiva.client import models
from inductiva.utils import format_utils


def _fetch_tasks_from_api(
    status: Optional[Union[str, models.TaskStatusCode]] = None,
    page=1,
    per_page=10,
) -> List[Dict]:
    """Get information about a user's tasks on the API.

    Tags can be filtered by a status. Results are paginated indexed from 1.
    """
    api_config = api.validate_api_key(inductiva.api_key)

    with ApiClient(api_config) as client:
        api_instance = TasksApi(client)

        query_params = {
            "page": page,
            "per_page": per_page,
        }

        if status is not None:
            query_params["status"] = models.TaskStatusCode(status)

        try:
            # Get User Tasks
            resp = api_instance.get_user_tasks(
                query_params=query_params,
                skip_deserialization=True,  # avoid deserializing to model,
                # leave as dict which we'll later serialize to our own
                # dataclasses
            ).response

            response_body = json.loads(resp.data.decode("utf-8"))

            return [{**task} for task in response_body]

        except ApiException as e:
            raise e


def _dataframe_of_tasks(tasks: Sequence["inductiva.tasks.Task"]) -> str:
    columns = [
        "ID", "Simulator", "Status", "Submitted", "Started", "Computation Time",
        "Resource Type"
    ]
    rows = []

    for task in tasks:
        info = task.get_info()
        simulator = task.get_simulator_name()
        status = task.get_status()

        computation_end_time = info.get("computation_end_time", None)

        execution_time = task.get_computation_time(fail_if_running=False)
        if execution_time is not None:
            if computation_end_time is None:
                if status in ["started", "submitted"]:
                    execution_time = f"*{execution_time}"
                else:
                    execution_time = "n/a"

        executer = info["executer"]
        if executer is None:
            resource_type = None
        else:
            resource_type = executer["vm_type"]
            if executer["n_mpi_hosts"] > 1:
                resource_type += f" x{executer['n_mpi_hosts']}"

        row = [
            task.id,
            simulator,
            status,
            format_utils.datetime_formatter(info.get("input_submit_time",
                                                     None)),
            format_utils.datetime_formatter(info.get("start_time", None)),
            execution_time,
            resource_type,
        ]
        rows.append(row)

    return pd.DataFrame(rows, columns=columns)


# pylint: disable=redefined-builtin
def list(last_n: int = 5,
         status: Optional[Union[str, models.TaskStatusCode]] = None) -> None:
    # pylint: disable=line-too-long
    """List the last N tasks of a user.

    This function returns a dataframe with information about the last N tasks
    (with respect to submission time) of a user, sorted by submission time with
    the most recent first.
    A status can be specified to filter to get only tasks with that status, in
    which case the last N tasks with that status will be listed.
    The number of tasks can be less than N if the aren't enough tasks that
    match the specified criteria.

    Example usage:
        # list the last 5 tasks that were successful
        listing = inductiva.tasks.list(5)
        # listing is a dataframe with 5 rows
        print(listing)
                           ID           Simulator         Status          Submitted            Started Computation Time        Resource Type
    i5bnjw9gznqom857o4ohos661 openfoam_foundation        started   15 Jan, 15:02:07   15 Jan, 15:02:26         *03m 18s     c2d-standard-112
    tdbs1viqiu0g83icnp048dtjj              reef3d        started   15 Jan, 14:41:16   15 Jan, 14:41:51         *23m 53s       c2-standard-60
    ca370nk18dz2enzvp64y8a1r8 openfoam_foundation         failed   15 Jan, 14:40:48   15 Jan, 14:43:16          03m 23s  c2d-standard-56 x 4
    rorv72bagv8s03qqoih5t9tba openfoam_foundation      submitted   15 Jan, 14:39:02                n/a              n/a                  n/a
    4yv6mcbyo8x6eewv2xdy2x8ws openfoam_foundation      submitted   15 Jan, 14:38:05                n/a              n/a                  n/a
    Args:
        last_n: The number of most recent tasks with respect to submission
            time to list. If filtering criteria (currently status is available)
            is specified, most recent N tasks that match that criteria will be
            listed. The actual number of tasks may be less if there
            aren't enough tasks available.
        status: The status of the tasks to list. If None, tasks with any status
            will be listed.
    """
    # pylint: enable=line-too-long
    status = models.TaskStatusCode(status) if status is not None else None
    tasks = get(last_n, status=status)
    tasks_df = _dataframe_of_tasks(tasks)

    return tasks_df


def get(
    last_n: int = 5,
    status: Optional[Union[str, models.TaskStatusCode]] = None
) -> List["inductiva.tasks.Task"]:
    """Get the last N tasks of a user.

    This function fetches info about the last N tasks (with respect to
    submission time) of a user to stdout, sorted by submission time with the
    most recent first.
    A status can be specified to filter to get only tasks with that status, in
    which case the last N tasks with that status will be listed.
    The number of tasks can be less than N if the aren't enough tasks that match
    the specified criteria.

    Similar to the inductiva.task.list() function, but instead of printing
    to stdout, returns a list of Task objects which can be used to perform
    further operations on those tasks (e.g., download outputs,
    kill unfinished tasks, etc.).

    Example usage:
        # get the last 5 tasks that haven't started yet and kill them
        tasks = inductiva.tasks.get(5, status="submitted")
        for task in tasks:
            task.kill()

    Args:
        last_n: The number of most recent tasks with respect to submission
            time to fetch. If filtering criteria (currently status is available)
            is specified, most recent N tasks that match that criteria will be
            listed. The actual number of tasks may be less if there
            aren't enough tasks available.
        status: The status of the tasks to get. If None, tasks with any status
            will be returned.

    Returns:
        List of Task objects.
    """
    if last_n < 1:
        raise ValueError("last_n must be >= 1")

    status = models.TaskStatusCode(status) if status is not None else None

    raw_tasks_info = _fetch_tasks_from_api(status, page=1, per_page=last_n)
    tasks = [
        inductiva.tasks.Task.from_api_info(info) for info in raw_tasks_info
    ]

    return tasks
