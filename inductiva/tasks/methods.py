"""Methods to interact with the tasks submitted to the API."""
import json
from typing import Dict, List, Optional, Union

from inductiva import tasks
from inductiva import client as api_client


def _fetch_tasks_from_api(
        status: Optional[str] = None,
        page=1,
        per_page=10,
        project_name: str = None) -> List[Dict]:
    """Get information about a user's tasks on the API.

    Tags can be filtered by a status. Results are paginated indexed from 1.
    """

    with api_client.methods.get_client() as client:
        api_instance = api_client.TasksApi(client)

        query_params = {
            "page": page,
            "per_page": per_page,
        }

        if project_name is not None:
            query_params["project"] = project_name

        if status is not None:
            query_params["status"] = api_client.models.TaskStatusCode(status)

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

        except api_client.ApiException as e:
            raise e



def get(
    last_n: int = None,
    status: Optional[str] = None,
    project: str = None,
) -> List[tasks.Task]:
    """Get tasks of a user.

    This function fetches info about the last N tasks of a user, sorted
    by submission time with the most recent first.
    A status can be specified to filter to get only tasks with that status, in
    which case the last N tasks with that status will be listed.
    The number of tasks can be less than N if the aren't enough tasks that match
    the specified criteria.

    Example usage:
        # get the last 5 tasks that haven't started yet and kill them
        tasks = inductiva.tasks.get(5, status="submitted")
        for task in tasks:
            task.kill()
    Args:
        last_n (int): The number of tasks with repect to the submission
            time to fectch. If `last_n` is `None` we fetch all tasks
            available mathching the remaining criteria.
        project (str): The project name from which to fetch. If None,
            fetches from all projects.
        status: Status of the tasks to get. If `None`, tasks with any
            status will be returned.
    Returns:
        List of dictionaries with information about the tasks.
    """

    if last_n and not(last_n >= 1):
        raise ValueError("if provided, last_n must be >= 1")

    all_tasks = _fetch_tasks_from_api(status,
                                                 page=1,
                                                 per_page=last_n,
                                                 project=project)
    page_counter = 1

    if last_n <= _MAX_TASKS_PER_PAGE:

    while tasks_fetched := _fetch_tasks_from_api(status,
                                                 page=page_counter,
                                                 per_page=500,
                                                 project=project):
        all_tasks.extend(tasks_fetched)
        page_counter += 1
        if len(all_tasks) >= last_n:
            all_tasks = all_tasks[:last_n]
            break

    tasks = [tasks.Task.from_api_info(info) for info in all_tasks]

    return tasks
