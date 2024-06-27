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
from inductiva.utils import format_utils


def compare(list_of_tasks: Iterable[Task],
            sort_by: str = "Computation",
            reverse: bool = False,
            verbose: int = 2) -> str:
    """Prints tasks sorted by a given metric and returns the ID of the first
    task.

    This method prints a table of tasks sorted by a given metric (computation by
    default). The sort order can be reversed. The ID of the first task in the
    sorted list is returned. This can represent the task with the lowest
    computation time, or the highest compression time, based on the sort_by and
    reverse parameters.

    Args:
        list_of_tasks: An iterable of tasks.
        sort_by: The metric to sort by. Allowed values: "Computation",
            "Compression", "Upload", "Total Cost".
        reverse: Whether to sort in reverse order.
        verbose: 0: no output, 1: only prints the table, 2: prints table and
        task info
    returns:
        The ID of the first task in the sorted dict. If the list is empty,
        returns None.
    """
    if not list_of_tasks:
        print("No tasks to compare.")
        return None

    final_table = defaultdict(list)

    str_to_index = {
        "ID": 0,
        "Computation": 1,
        "Compression": 2,
        "Upload": 3,
        "Total Cost": 4,
    }

    if sort_by not in str_to_index:
        raise ValueError(f"Invalid sort_by value: {sort_by}.\n"
                         f"Allowed values: {list(str_to_index.keys())}")

    for task in list_of_tasks:
        if not task.info.status == models.TaskStatusCode.SUCCESS:
            continue
        final_table["ID"].append(task.id)

        final_table["Computation"].append(
            task.info.time_metrics.computation_seconds.value)
        final_table["Compression"].append(
            task.info.time_metrics.output_compression_seconds.value)
        final_table["Upload"].append(
            task.info.time_metrics.output_upload_seconds.value)
        final_table["Total cost"].append(1.0)

    combined_data = list(zip(*final_table.values()))

    sorted_data = sorted(combined_data,
                         key=lambda x: x[str_to_index[sort_by]],
                         reverse=reverse)

    # Unpack the sorted data back into the defaultdict
    sorted_data_dict = defaultdict(list)
    for entry in sorted_data:
        sorted_data_dict["ID"].append(entry[0])
        sorted_data_dict["Computation"].append(entry[1])
        sorted_data_dict["Compression"].append(entry[2])
        sorted_data_dict["Upload"].append(entry[3])
        sorted_data_dict["Total Cost"].append(entry[4])

    emph_formatter = format_utils.get_ansi_formatter()

    header_formatters = [
        lambda x: emph_formatter(x.upper(), format_utils.Emphasis.BOLD)
    ]

    list_str_sorted = format_utils.get_tabular_str(
        sorted_data_dict, header_formatters=header_formatters)

    first_task_id = sorted_data_dict["ID"][0]

    if verbose > 0:
        print(f"Tasks sorted by {sort_by}:")
        print(list_str_sorted)
        if verbose > 1:
            print(f"The first task id is -> {first_task_id}")
            print(inductiva.tasks.Task(first_task_id).get_info())

    return first_task_id


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
            resource_type = (f"{task.info.executer.host_type} "
                             f"{task.info.executer.vm_type}")
            if task.info.executer.n_mpi_hosts > 1:
                resource_type += f" x{task.info.executer.n_mpi_hosts}"

        table["ID"].append(task.id)
        table["Simulator"].append(task.get_simulator_name())
        table["Status"].append(task.info.status)
        table["Submitted"].append(task.info.input_submit_time)
        table["Started"].append(task.info.start_time)
        table["Computation Time"].append(execution_time)
        table["Resource Type"].append(resource_type)

    return table


def _fetch_tasks_from_api(
        status: Optional[Union[str, models.TaskStatusCode]] = None,
        page=1,
        per_page=10,
        project: Union[str, "inductiva.projects.Project"] = None) -> List[Dict]:
    """Get information about a user's tasks on the API.

    Tags can be filtered by a status. Results are paginated indexed from 1.
    """
    api_config = api.get_api_config()

    with ApiClient(api_config) as client:
        api_instance = TasksApi(client)

        query_params = {
            "page": page,
            "per_page": per_page,
        }

        if project is not None:
            if isinstance(project, inductiva.projects.Project):
                project = project.name
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


def get(
    last_n: int = 5,
    status: Optional[Union[str, models.TaskStatusCode]] = None,
    project: Union[str, "inductiva.projects.Project"] = None
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
        project: The project from which to fetch. If None, fetches from all
            projects.

    Returns:
        List of Task objects.
    """
    if last_n < 1:
        raise ValueError("last_n must be >= 1")

    status = models.TaskStatusCode(status) if status is not None else None

    raw_tasks_info = _fetch_tasks_from_api(status,
                                           page=1,
                                           per_page=last_n,
                                           project=project)
    tasks = [
        inductiva.tasks.Task.from_api_info(info) for info in raw_tasks_info
    ]

    return tasks


def get_all(
    status: Optional[Union[str, models.TaskStatusCode]] = None,
    project: Union[str, "inductiva.projects.Project"] = None,
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
    status = models.TaskStatusCode(status) if status is not None else None

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
