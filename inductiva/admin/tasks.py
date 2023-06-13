"""Admin functions to query tasks."""
from typing import Dict, List, Optional

import inductiva
from inductiva import api
from inductiva.client import ApiClient, ApiException
from inductiva.client.apis.tags.admin_api import AdminApi
from inductiva.client.model.task_status_code import TaskStatusCode


def get_tasks_info(statuses: Optional[List[str]] = None,
                   usernames: Optional[List[str]] = None,
                   methods: Optional[List[str]] = None,
                   page=1,
                   per_page=10) -> List[Dict]:
    """Get information about tasks all tasks on the API.

    Tasks can be filtered by status, username, and method_name,
    by providing a list of each as function arguments. Results are
    paginated indexed from 1.
    """
    api_config = api.validate_api_key(inductiva.api_key)

    with ApiClient(api_config) as client:
        api_instance = AdminApi(client)

        query_params = {
            "page": page,
            "per_page": per_page,
        }
        if statuses:
            query_params["status"] = [
                TaskStatusCode(status) for status in statuses
            ]

        if usernames:
            query_params["username"] = usernames

        if methods:
            query_params["method"] = methods

        try:
            # Get User Tasks
            api_response = api_instance.get_tasks(query_params=query_params,)

            return [{**task} for task in api_response.body]

        except ApiException as e:
            raise e
