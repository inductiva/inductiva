"""Methods to interact with the tasks submitted to the API."""
import json
from typing import Dict, List, Optional, Union, Sequence

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


def _list_of_tasks_to_str(tasks: Sequence["inductiva.tasks.Task"]) -> str:
    columns = [
        "ID", "Scenario", "Simulator", "Status", "Submitted", "Started",
        "Duration", "VM Type"
    ]
    rows = []

    for task in tasks:
        info = task.get_info()
        # e.g., get "openfoam" from "fvm.openfoam.run_simulation"
        simulator = task.get_simulator_name()
        scenario = task.get_scenario_name()
        if scenario:
            scenario = scenario.replace("_", " ")

        end_time = info.get("end_time", None)
        execution_time = task.get_execution_time(fail_if_running=False)

        if execution_time is not None:
            execution_time = format_utils.seconds_formatter(execution_time)
            if end_time is None:
                execution_time = f"*{execution_time}"

        row = [
            task.id,
            scenario,
            simulator,
            task.get_status(),
            info.get("input_submit_time", None),
            info.get("start_time", None),
            execution_time,
            task.get_machine_type(),
        ]
        rows.append(row)
    formatters = {
        "Submitted": format_utils.datetime_formatter,
        "Started": format_utils.datetime_formatter,
    }

    override_col_space = {
        "Scenario": 22,
        "Submitted": 20,
        "Started": 20,
        "Status": 20,
        "VM Type": 18,
    }

    return format_utils.get_tabular_str(
        rows,
        columns,
        default_col_space=15,
        override_col_space=override_col_space,
        formatters=formatters,
    )


# pylint: disable=redefined-builtin
def list(last_n: int = 5,
         status: Optional[Union[str, models.TaskStatusCode]] = None) -> None:
    # pylint: disable=line-too-long
    """List the last N tasks of a user.

    This function lists info about the last N tasks (with respect to submission
    time) of a user to stdout, sorted by submission time with the
    most recent first.
    A status can be specified to filter to get only tasks with that status, in
    which case the last N tasks with that status will be listed.
    The number of tasks can be less than N if the aren't enough tasks that
    match the specified criteria.

    Example usage:
        # list the last 5 tasks that were successful
        inductiva.tasks.list(5, status="success")
                         ID               Scenario       Simulator               Status            Submitted              Started        Duration            VM Type
        1695309310050950437                    n/a    dualsphysics              started     21 Sep, 16:15:10     21 Sep, 16:15:10       *0h 2m 5s      c2-standard-4
        1695307352995292367      protein solvation         gromacs              started     21 Sep, 15:42:33     21 Sep, 15:42:35     *0h 34m 41s      c2-standard-4
        1695306097851999678      protein solvation         gromacs              success     21 Sep, 15:21:38     21 Sep, 15:21:40      0h 15m 43s      c2-standard-4
        1695294928922012701            wind tunnel        openfoam              success     20 Sep, 12:15:31     20 Sep, 12:15:33       0h 0m 50s      c2-standard-4
        1695294363618736570            wind tunnel        openfoam              success     20 Sep, 12:06:06     20 Sep, 12:06:06       0h 0m 50s      c2-standard-4

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
    tasks_str = _list_of_tasks_to_str(tasks)
    print(tasks_str)


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
