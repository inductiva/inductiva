"""Methods to interact with the user storage resources."""
import inductiva
from inductiva.client.apis.tags import instance_api
from inductiva.utils import format_utils
from typing import Literal


def get_space_used():
    """Returns the occupied storage size in GB."""
    try:
        api = instance_api.InstanceApi(inductiva.api.get_client())
        response = api.get_storage_size()
        storage_used = round(float(response.body), 3)
        print(f"Total storage used: {storage_used} GB")
        return storage_used
    except inductiva.client.ApiException as api_exception:
        raise api_exception


def list_contents(max_results: int = 3,
                  sort_by: Literal["size", "creation_time"]= "size",
                  order: Literal["asc", "desc"] = "desc"):
    """List and display the contents of the user storage.

    Args:
        max_results (int): The maximum number of results to return.
        sort_by (str): The field to sort the contents by.
        order (str): Whether to sort the contents in ascending or
        descending order.
    Returns:
        list of dict: A list of dictionaries containing information about 
        the size, the name and the creation time of each content.

    This function prints a table with the storage content information: 
        Name            Size            Creation Time
        1234            5.68 MiB        29 Sep, 14:12:00
        12345           374.85 KiB      29 Sep, 14:13:10
        1234567         97.59 KiB       29 Sep, 14:13:24
        123             0 B             29 Sep, 14:13:29
        123456          0 B             29 Sep, 14:13:18
    """
    try:
        api = instance_api.InstanceApi(inductiva.api.get_client())
        contents = api.list_storage_contents({
            "max_results": max_results,
            "sort_by": sort_by,
            "order": order
        }).body
        all_contents = []
        for content_name, info in contents.items():
            size = info["size_bytes"]
            creation_time = info["creation_time"]
            all_contents.append({
                "content_name": content_name,
                "size": round(float(size), 3),
                "creation_time": creation_time
            })
        print(_print_contents_table(all_contents))
        return all_contents
    except inductiva.client.ApiException as api_exception:
        raise api_exception


def _print_contents_table(contents):
    columns = ["Name", "Size", "Creation Time"]
    rows = []

    for content in contents:
        row = [
            content['content_name'], content['size'], content['creation_time']
        ]
        rows.append(row)
    override_col_space = {
        "Name": 22,
        "Size": 20,
        "Creation Time": 20,
    }
    formatters = {
        "Creation Time": format_utils.datetime_formatter,
        "Size": format_utils.bytes_formatter
    }

    return format_utils.get_tabular_str(
        rows,
        columns,
        default_col_space=15,
        override_col_space=override_col_space,
        formatters=formatters,
    )

def delete_task_storage(task_id: str):
    """Deletes the storage associated with a task.
        Args: 
            task_id (str): The task id."""
    try:
        api = instance_api.InstanceApi(inductiva.api.get_client())
        response = api.delete_task_storage({"task_id": task_id})
        print(response.body)
    except inductiva.client.ApiException as api_exception:
        raise api_exception 

