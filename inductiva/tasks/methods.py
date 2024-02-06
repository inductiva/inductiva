"""Methods to interact with the tasks submitted to the API."""
from collections import defaultdict
import json
from typing import Any, Dict, Iterable, List, Mapping, Optional, Union

import inductiva
from inductiva import api
from inductiva.client import ApiClient, ApiException
from inductiva.client.apis.tags.tasks_api import TasksApi
from inductiva.client import models
from inductiva.tasks.task import Task


def to_dict(list_of_tasks: Iterable[Task]) -> Mapping[str, List[Any]]:
    """
    Converts an Iterable of tasks to a dictionary with all the
    relevant information for all the tasks.
        Args:
            list_of_tasks: An Iterable of tasks.
        Returns:
            A dictionary with all the relevant information for
            all the tasks. Example: { "ID": [1, 2, 3], 
            "Simulator": ["reef3d", "reef3d", "reef3d"], ... }
    """

    table = defaultdict(list)

    for task in list_of_tasks:
        info = task.get_info()
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
            n_mpi_hosts = executer["n_mpi_hosts"]
            if n_mpi_hosts > 1:
                resource_type += f" x{n_mpi_hosts}"
        table["ID"].append(task.id)
        table["Simulator"].append(task.get_simulator_name())
        table["Status"].append(status)
        table["Submitted"].append(info.get("input_submit_time", None))
        table["Started"].append(info.get("start_time", None))
        table["Computation Time"].append(execution_time)
        table["Resource Type"].append(resource_type)
    return table


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


def get_all(
        status: Optional[Union[str,
                               models.TaskStatusCode]] = None) -> List[Dict]:
    """Get all tasks of a user.

    This function fetches all tasks of a user, sorted by submission
    time with the most recent first. If status is specified, only
    tasks with that status will be fetched.
    Args:
        status: The status of the tasks to get. If None, tasks with any status
            will be returned.
    Returns:
        List of dictionaries with information about the tasks.
    """
    status = models.TaskStatusCode(status) if status is not None else None

    all_tasks = []
    page_counter = 1

    while tasks_fetched := _fetch_tasks_from_api(status,
                                                 page=page_counter,
                                                 per_page=50):
        all_tasks.extend(tasks_fetched)
        page_counter += 1
    return all_tasks
