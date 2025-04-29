"""Methods to interact with the tasks submitted to the API."""
from collections import defaultdict
import json
from typing import Any, Dict, Iterable, List, Mapping, Optional

import inductiva
from inductiva.client import models
from inductiva.tasks.task import Task
from inductiva.utils import format_utils
from inductiva.client import ApiException
from inductiva.client.apis.tags.tasks_api import TasksApi


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
    column_names = [
        "ID", "Simulator", "Status", "Submitted", "Started", "Computation Time",
        "Resource Type"
    ]
    table = defaultdict(list, {key: [] for key in column_names})

    for task in list_of_tasks:
        execution_time = task.get_computation_time(cached=True)

        if execution_time is not None:
            execution_time = format_utils.seconds_formatter(execution_time)
            if task.info.computation_end_time is None:
                if task.info.status in ["started", "submitted"]:
                    execution_time = f"*{execution_time}"
                else:
                    execution_time = "n/a"

        if task.info.executer is None:
            resource_type = None
        else:
            if task.info.executer.vm_type == "n/a":
                vm_type = task.info.executer.vm_name
            else:
                vm_type = task.info.executer.vm_type
            resource_type = (f"{task.info.executer.host_type} "
                             f"{vm_type}")
            if task.info.executer.n_mpi_hosts > 1:
                resource_type += f" x{task.info.executer.n_mpi_hosts}"

        table["ID"].append(task.id)
        table["Simulator"].append(task.get_simulator_name())
        table["Status"].append(task.info.status_alias)
        table["Submitted"].append(task.info.input_submit_time)
        table["Started"].append(task.info.start_time)
        table["Computation Time"].append(execution_time)
        table["Resource Type"].append(resource_type)

    return table


def _fetch_tasks_from_api(status: Optional[str] = None,
                          page=1,
                          per_page=10,
                          project: Optional[str] = None) -> List[Dict]:
    """Get information about a user's tasks on the API.

    Tags can be filtered by a status. Results are paginated indexed from 1.
    """

    with inductiva.api.methods.get_client() as client:
        api_instance = TasksApi(client)

        query_params = {
            "page": page,
            "per_page": per_page,
        }

        if project is not None:
            query_params["project"] = project

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


def get(last_n: int = 5,
        status: Optional[str] = None,
        project: Optional[str] = None) -> List["inductiva.tasks.Task"]:
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
        project: The project from which to fetch. If None, fetches from all
            projects.

    Returns:
        List of Task objects.
    """
    if last_n < 1:
        raise ValueError("last_n must be >= 1")

    raw_tasks_info = _fetch_tasks_from_api(status,
                                           page=1,
                                           per_page=last_n,
                                           project=project)
    tasks = [
        inductiva.tasks.Task.from_api_info(info) for info in raw_tasks_info
    ]

    return tasks


def get_tasks(last_n: int = 10,
              project: Optional[str] = None,
              status: Optional[str] = None):
    """Get the last N submitted tasks.

        Get the last N submitted tasks, eventually filtered by status.
        By default, only the last 10 submitted tasks are returned,
        irrespectively of their status.

        Args:
            last_n (int): The number of tasks with repect to the submission
                time to fectch. If `last_n<=0` we fetch all tasks submitted
                to the project.
            status: Status of the tasks to get. If `None`, tasks with any
                status will be returned.
        """
    if last_n <= 0:
        return inductiva.tasks.get_all(status=status, project=project)
    return inductiva.tasks.get(last_n=last_n, status=status, project=project)


def get_all(
    status: Optional[str] = None,
    project: Optional[str] = None,
) -> List["inductiva.tasks.Task"]:
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
    all_tasks = []
    page_counter = 1
    while tasks_fetched := _fetch_tasks_from_api(status,
                                                 page=page_counter,
                                                 per_page=500,
                                                 project=project):
        all_tasks.extend(tasks_fetched)
        page_counter += 1

    tasks = [inductiva.tasks.Task.from_api_info(info) for info in all_tasks]

    return tasks
