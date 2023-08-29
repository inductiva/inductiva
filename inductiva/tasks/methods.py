"""Methods to interact with the tasks submitted to the API."""
import json
from typing import Dict, List, Optional, Union

import inductiva
from inductiva import api
from inductiva.client import ApiClient, ApiException
from inductiva.client.apis.tags.tasks_api import TasksApi
from inductiva.client import models


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


# pylint: disable=redefined-builtin
def list(num_tasks,
         status: Optional[Union[str, models.TaskStatusCode]] = None) -> None:
    # pylint: disable=line-too-long
    """List the last N tasks of a user.

    This function will fetch info on the last N tasks of a user and print
    them to stdout formatted as a table. A status can be specified to filter
    to get only tasks with that status. Tasks are sorted by submission time
    with the most recent first.

    Example usage:
        # list the last 5 tasks that were successful
        inductiva.tasks.list(5, status="success")
        ID                   Simulator    Status       Submitted            Started              Duration     VM Type
        1691150776862178362  openfoam     success      04 Aug, 12:06:17     04 Aug, 12:06:18     0h 1m 53s    c2-standard-8
        1691149904961476240  openfoam     success      04 Aug, 11:51:46     04 Aug, 11:51:46     0h 1m 28s    c2-standard-8
        1691081881158823776  openfoam     success      03 Aug, 16:58:02     03 Aug, 16:58:02     0h 1m 20s    n2-standard-32
        1691081409916619414  openfoam     success      03 Aug, 16:50:11     03 Aug, 16:50:11     0h 1m 20s    n2-standard-32
        1691080520213617518  openfoam     success      03 Aug, 16:35:21     03 Aug, 16:35:21     0h 1m 23s    n2-standard-32

    Args:
        num_tasks: The number of tasks to list.
        status: The status of the tasks to list. If None, tasks with any status
            will be listed.
    """
    # pylint: enable=line-too-long
    status = models.TaskStatusCode(status) if status is not None else None

    id_header = "ID"
    simulator_header = "Simulator"
    status_header = "Status"
    submitted_header = "Submitted"
    started_header = "Started"
    duration_header = "Duration"
    vm_type_header = "VM Type"

    print(f"  {id_header:<20} {simulator_header:<12} {status_header:<12}"
          f" {submitted_header:<20} {started_header:<20} {duration_header:<12}"
          f" {vm_type_header:<12}")

    tasks = get(num_tasks, status=status)
    for task in tasks:
        print(task)


def get(
    num_tasks,
    status: Optional[Union[str, models.TaskStatusCode]] = None
) -> List["inductiva.tasks.Task"]:
    """Get the last N tasks of a user.

    This function will fetch info on the last N tasks of a user and return
    a list of Task objects. A status can be specified to filter to get only
    tasks with that status. Tasks are sorted by submission time with the
    most recent first.

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
        num_tasks: The number of tasks to get.
        status: The status of the tasks to get. If None, tasks with any status
            will be returned.

    Returns:
        List of Task objects.
    """
    status = models.TaskStatusCode(status) if status is not None else None

    raw_tasks_info = _fetch_tasks_from_api(status, page=1, per_page=num_tasks)
    tasks = [
        inductiva.tasks.Task.from_api_info(info) for info in raw_tasks_info
    ]

    return tasks
