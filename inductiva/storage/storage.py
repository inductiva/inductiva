"""Methods to interact with the user storage resources."""
import inductiva
from inductiva.client.apis.tags import storage_api
from inductiva.utils import format_utils
from typing import Literal


def get_space_used():
    """Returns the occupied storage size in GB."""
    try:
        api = storage_api.StorageApi(inductiva.api.get_client())
        response = api.get_storage_size()
        storage_used = round(float(response.body), 3)
        print(f"Total storage used: {storage_used} GB")
        return storage_used
    except inductiva.client.ApiException as api_exception:
        raise api_exception


def listdir(path="/",
            max_results: int = 10,
            order_by: Literal["size", "creation_time"] = "size",
            sort_order: Literal["asc", "desc"] = "desc"):
    """List and display the contents of the user's storage.
    Args:
        path (str): Storage directory to list. Default is root.
        max_results (int): The maximum number of results to return.
        order_by (str): The field to sort the contents by.
        sort_order (str): Whether to sort the contents in ascending or
        descending order.
    Returns:
        list of dict: A list of dictionaries containing information about
        the size, the name and the creation time of each content that can
        easily be converted to a dataframe.

    This function prints a table with the storage content information:
        Name            Size            Creation Time
        1234            5.68 MiB        29 Sep, 14:12:00
        12345           374.85 KiB      29 Sep, 14:13:10
        1234567         97.59 KiB       29 Sep, 14:13:24
        123             0 B             29 Sep, 14:13:29
        123456          0 B             29 Sep, 14:13:18
    You can use this information to delete the contents you don't need
    anymore and further inspect task outputs and logs using the Task
    class.
    """
    try:
        api = storage_api.StorageApi(inductiva.api.get_client())
        contents = api.list_storage_contents({
            "dir_name": path,
            "max_results": max_results,
            "sort_by": order_by,
            "order": sort_order
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
            content["content_name"], content["size"], content["creation_time"]
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


def rmdir(path: str):
    """Deletes a directory inside the user's storage in Inductiva API.
    Currently, directories are named by the task id of the task that
    created them, so this function can be used to delete the directory
    of a task.
    Args:
        path (str): The path to the directory to delete.
    """
    try:
        api = storage_api.StorageApi(inductiva.api.get_client())
        api.delete_directory({"dir_name": path})
        stripped_path = path.rstrip("/")
        print(f"Directory deleted: {stripped_path}")
    except inductiva.client.ApiException as api_exception:
        raise api_exception
