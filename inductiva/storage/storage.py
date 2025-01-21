"""Methods to interact with the user storage resources."""
import time
import logging
import os
import tqdm
from typing import List, Literal, Optional, Tuple
import urllib

import inductiva
from inductiva import constants
from inductiva.api import methods
from inductiva.client import exceptions, models
from inductiva.client.apis.tags import storage_api
from inductiva.utils import format_utils


def _print_storage_size_and_cost() -> int:
    """ Print the storage size and cost.

    return the storage size in bytes.
    """
    api = storage_api.StorageApi(inductiva.api.get_client())
    storage_total_size_bytes = api.get_storage_size().body
    estimated_storage_cost = api.get_storage_monthly_cost(
    ).body["estimated_monthly_cost"]

    estimated_storage_cost = format_utils.currency_formatter(
        estimated_storage_cost)
    storage_total_size = format_utils.bytes_formatter(storage_total_size_bytes)

    print("Total storage size used:")
    print(f"\tVolume: {storage_total_size}")
    print(f"\tCost: {estimated_storage_cost}/month")
    print("")

    return storage_total_size_bytes


def get_space_used():
    """Returns the occupied storage size in GB."""
    storage_total_size_bytes = _print_storage_size_and_cost()
    #Return float instead of a string
    storage_used = round(float(storage_total_size_bytes) / (1024**3), 3)
    return storage_used


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

    api = storage_api.StorageApi(inductiva.api.get_client())

    # This is valid for single files and directories
    if len(path.split(os.sep)) < 2:
        path += os.sep

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

    _print_storage_size_and_cost()

    print(f"Listed {len(all_contents)} folder(s). Ordered by {order_by}.\n"
          "Use --max-results/-m to control the number of results displayed.")
    return all_contents


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


def upload_from_url(
    url: str,
    remote_dir: str,
    file_name: Optional[str] = None,
):
    """
    Upload a file from a given URL to a specified remote directory.

    If no file name is provided, the function extracts the name from the URL.

    Args:
        url (str): The URL of the file to be uploaded.
        remote_dir (str): The path to the remote directory where the file will
            be stored.
        file_name (str, optional): The name to save the uploaded file as.
            If not specified, the name will be extracted from the URL.
    """
    api_instance = storage_api.StorageApi(inductiva.api.get_client())

    # append filename extracted from url to remote_path
    if not file_name:
        parsed_url = urllib.parse.urlparse(url)
        file_name = urllib.parse.unquote(parsed_url.path.split("/")[-1])

    remote_path = os.path.join(remote_dir, file_name)

    api_instance.upload_from_url(query_params={
        "url": url,
        "path": remote_path,
        "unzip": "f",
    },)
    logging.info("File is being uploaded...")
    logging.info("You can use 'inductiva storage ls' to check the status.")


def upload(
    local_path: str,
    remote_dir: str,
):
    """
    Upload a local file or directory to a specified remote directory.

    Args:
        local_path (str): The path to the local file or directory to be
            uploaded.
        remote_dir (str, optional): The remote directory where the file will
            be uploaded.
    """
    is_dir = os.path.isdir(local_path)

    if is_dir:
        local_dir = os.path.join(local_path, "")
        file_paths, total_size = _list_files(local_path)
        remote_file_paths = [
            os.path.join(remote_dir, file_path) for file_path in file_paths
        ]
    else:
        local_dir = os.path.dirname(local_path)
        filename = os.path.basename(local_path)
        remote_file_paths = [os.path.join(remote_dir, filename)]
        total_size = os.path.getsize(local_path)

    if os.path.join(remote_dir, constants.TASK_OUTPUT_ZIP) in remote_file_paths:
        raise ValueError(f"Invalid file name: '{constants.TASK_OUTPUT_ZIP}.'")

    logging.info("Uploading input...")

    api_instance = storage_api.StorageApi(inductiva.api.get_client())

    api_response = methods.get_upload_url(
        api_instance.get_upload_url,
        query_params={"paths": remote_file_paths},
    )

    with tqdm.tqdm(total=total_size,
                   unit="B",
                   unit_scale=True,
                   unit_divisor=1000) as progress_bar:

        for response in api_response:
            method = response["method"]
            url = response["url"]
            remote_file_path = response["file_path"]

            file_path = remote_file_path.removeprefix(f"{remote_dir}/")
            local_file_path = os.path.join(local_dir, file_path)

            try:
                methods.upload_file(api_instance, local_file_path, method, url,
                                    progress_bar)

                methods.notify_upload_complete(
                    api_instance.notify_upload_file,
                    query_params={
                        "path": remote_file_path,
                    },
                )
            except exceptions.ApiException as e:
                raise e

    logging.info("Input uploaded successfully.")


def _list_files(root_path: str) -> Tuple[List[str], int]:
    """
    Lists all files within a directory and returns their total size.

    Args:
        root_path: The root directory to list files from.
    """
    file_paths = []
    total_size = 0

    for dirpath, _, filenames in os.walk(root_path):
        for filename in filenames:
            path = os.path.relpath(os.path.join(dirpath, filename), root_path)
            file_paths.append(path)

            full_path = os.path.join(dirpath, filename)
            total_size += os.path.getsize(full_path)

    return file_paths, total_size


def remove_workspace(remote_dir) -> bool:
    """
    Removes path from a remote directory.

    Parameters:
    - remote_dir (str): The path to the remote directory.
    """
    api = storage_api.StorageApi(inductiva.api.get_client())

    logging.info("Removing workspace file(s)...")

    # Since we don't allow root files in workspaces it must be a directory
    # otherwise path validation in the backend will give error
    if "/" not in remote_dir:
        remote_dir = remote_dir + "/"
    api.delete_file(query_params={"path": remote_dir},)
    logging.info("Workspace file(s) removed successfully.")


class StorageOperation():
    """Represents a storage operation running remotely via Inductiva API."""

    def __init__(self, api, id_):
        self._api = api
        self.id = id_

    def _update_from_api_response(self, response):
        self._name = response["name"]
        self._status = response["status"]
        self._attributes = response["attributes"]
        self._start_time = response["start_time"]
        self._end_time = response["end_time"]
        self._error_message = response["error_message"]

    @classmethod
    def from_api_response(cls, api, response):

        op = cls(api, response["id"])

        op._update_from_api_response(response,)
        return op

    def _refresh(self):
        resp = self._api.get_operation(path_params={
            "operation_id": self.id
        }).body

        self._update_from_api_response(resp)

    def wait(self, poll_s: int = 2):
        """Wait for the operation to complete.

        Args:
            poll_s: Time in seconds between calls to the API to update
                the status of the operation.
        """

        while self._status == models.OperationStatus.RUNNING:
            self._refresh()
            time.sleep(poll_s)

        if self._status == models.OperationStatus.FAILED:
            logging.error("Operation failed: %s", self._error_message)
        else:
            logging.info("Operation completed successfully.")

        return self._status


def export(path: str, dest_url: str) -> StorageOperation:
    """Export files from the API's storage to a remote storage location.

    Args:
        path: Path in the API's remote storage to the files.
        dest_url: URL to upload the output files.

    Returns:
        Instance of StorageOperation. Call `wait` on the resulting object
        to block until the operation finishes.
    """
    api = storage_api.StorageApi(inductiva.api.get_client())
    resp = api.export_files(body={
        "path": path,
        "dest_url": dest_url,
    })

    logging.info("Started export operation ...")
    operation = StorageOperation.from_api_response(api, resp.body)

    return operation
