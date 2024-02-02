"""Methods to interact with the user storage resources."""
from absl import logging

import inductiva
from inductiva.client import exceptions
from inductiva.client.apis.tags import storage_api
from inductiva.utils import format_utils
from typing import Literal


def get_space_used():
    """Returns the occupied storage size in GB."""
    try:
        api = storage_api.StorageApi(inductiva.api.get_client())
        response = api.get_storage_size()
        storage_used = round(float(response.body), 3)
        print(f"Total user's remote storage in use: {storage_used} GB")
        return storage_used
    except inductiva.client.ApiException as api_exception:
        raise api_exception


def listdir(path="/",
            max_results: int = 10,
            order_by: Literal["size", "creation_time"] = "creation_time",
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
            "path": path,
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

    formatters = {
        "Creation Time": [format_utils.datetime_formatter],
        "Size": [format_utils.bytes_formatter],
    }

    emph_formatter = format_utils.get_ansi_formatter()
    header_formatters = [
        lambda x: emph_formatter(x.upper(), format_utils.Emphasis.BOLD)
    ]

    return format_utils.get_tabular_str(
        rows,
        columns,
        formatters=formatters,
        header_formatters=header_formatters,
    )


def rmdir(path: str, /, confirm: bool = False):
    """Delete paths inside the user's remote storage in the Inductiva platform.

    This function removes the given `path` from the user's remote storage
    in the Inductiva platform. The given path can either be the storage name of
    a task to remove all content pertaining to that task, or `*` to remove all
    content. The function requires an explicit `confirm=True` keyword argument
    to confirm that the path is, indeed, to be removed.

    E.g.:
        #Works - remote path "task_id" exists and user is explicitly
        # sets confirm=True:
        >>> inductiva.storage.rmdir("task_id", confirm=True)
        
        #Works - The user explicitly sets confirm=True, confirming that
        # all contents in his remote storage space are to be removed:
        >>> inductiva.storage.rmdir("*", confirm=True)
        # Fails -  remote path "task_id" exists but the user fails to
        # confirm his intent. The call will fail win a runtime exception:
        >>> inductiva.storage.rmdir("task_id")

    Args:
        path (str): Path relative to the root of the bucket to delete.
            If path="*" then all the contents in the user's remote storage space
            will be removed.
        confirm (bool): Confirmation flag to explicitly ensure the intent
            to remove the given path. The call will fail with a RuntimeException
            if the user does not set `confirm=True`.
    """
    if not confirm:
        if path == "*":
            path = "all contents"
        raise RuntimeError("Please set `confirm=True` to confirm you want to "
                           f"delete {path} from the user's remote storage.")

    if path == "*":
        logging.info("Removing everything from user's remote storage.")
    else:
        logging.info("Removing %s in the user's remote storage.", path)

    try:
        api = storage_api.StorageApi(inductiva.api.get_client())
        api.delete_path({"path": path})
        logging.info("Successfully removed remote path '%s'.", path)
    except exceptions.ApiException as api_exception:
        if api_exception.status == 404:
            raise exceptions.ApiValueError(
                f"Unable to remove path '{path}'. Path does not exist in "
                "user's remote storage.")
        elif api_exception.status == 500:
            raise exceptions.ApiRuntimeError(
                f"Failed to remove remote path '{path}'.")
        raise api_exception
