"""Methods to interact with the user storage resources."""
import json
import logging
import os
import tempfile
import zipfile
import tqdm
import urllib3
from typing import List, Literal
from urllib.parse import unquote, urlparse

import inductiva
from inductiva import constants
from inductiva.client import exceptions, ApiException
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
    try:
        storage_total_size_bytes = _print_storage_size_and_cost()
        #Return float instead of a string
        storage_used = round(float(storage_total_size_bytes) / (1024**3), 3)
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
    remote_path: str = None,
    unzip: bool = False,
):
    """
    Upload a file from a URL to the user workspace.

    Args:
        url (str): The URL of the file to upload.
        remote_dir (str): The remote directory to upload the file to. 
        remote_path (str, optional): The path to save the file as. If not
            provided, the path will be extracted from the URL.
        unzip (bool, optional): Whether to unzip the file after uploading.
            Default is False.
    """
    api_instance = storage_api.StorageApi(inductiva.api.get_client())

    if remote_path is None:
        parsed_url = urlparse(url)
        remote_path = unquote(parsed_url.path.split("/")[-1])
    contents = api_instance.upload_from_url(
        query_params={
            "url": url,
            "file_path": remote_path,
            "unzip": "t" if unzip else "f",
        },
        path_params={
            "folder_name": remote_dir,
        },
    )
    if contents.response.status == 200:
        logging.info("File is being uploaded...")
    else:
        logging.info("File upload failed.")


def upload(
    local_path: str,
    remote_dir: str,
    remote_path: str = None,
):
    """
    Upload a local file or directory to the user workspace.

    Args:
        local_path (str): The path to the local file or directory to be
            uploaded.
        remote_dir (str, optional): The remote directory where the file will
            be uploaded. Defaults to "default".
        remote_path (str, optional): The remote path to save the file as.
    """
    api_instance = storage_api.StorageApi(inductiva.api.get_client())

    prefix_path = os.path.dirname(
        remote_path.lstrip("/")) if remote_path else ""

    is_dir = os.path.isdir(local_path)
    if is_dir:
        input_files = _list_files(local_path, prefix_path)
    else:
        file_name = os.path.basename(local_path)
        input_files = [os.path.join(prefix_path, file_name)]

    logging.info("Uploading input...")

    api_response = api_instance.get_upload_url(
        query_params={
            "file_paths": input_files,
        },
        path_params={
            "folder_name": remote_dir,
        },
    )

    for item in api_response.body:
        method = item["method"]
        url = item["url"]
        file_server_available = bool(item["file_server_available"])
        file_path = item["file_path"]
        file_path_local = file_path.removeprefix(prefix_path)
        if is_dir:
            file_path_local = os.path.join(local_path, file_path_local)
        else:
            file_path_local = os.path.join(os.path.dirname(local_path),
                                           file_path_local)

        headers = {"Content-Type": "application/octet-stream"}

        if not file_server_available:
            headers[
                "X-API-Key"] = api_instance.api_client.configuration.api_key[
                    "APIKeyHeader"]

        logging.debug("Upload URL: %s", url)

        file_size = os.path.getsize(file_path_local)

        with open(file_path_local, "rb") as zip_fp, tqdm.tqdm(
                total=file_size,
                unit="B",
                unit_scale=True,
                unit_divisor=1000,
        ) as progress_bar:
            # Wrap the file object so that the progress bar is updated
            # every time a chunk is read.
            wrapped_file = tqdm.utils.CallbackIOWrapper(
                progress_bar.update,
                zip_fp,
                "read",
            )
            # Use the pool_manager from the API client to send the request
            pool_manager: urllib3.PoolManager = (
                api_instance.api_client.rest_client.pool_manager)

            resp = pool_manager.request(method,
                                        url,
                                        body=wrapped_file,
                                        headers=headers)
            if resp.status != 200:
                raise ApiException(status=resp.status, reason=resp.reason)

            # Notify the API that the file upload is complete
            api_instance.notify_upload_file(
                query_params={
                    "file_path": file_path,
                    "unzip": "f",
                },
                path_params={
                    "folder_name": remote_dir,
                },
            )


def _zip_folder(source_path):
    """
    Zips a file or a folder and saves it to a temporary folder.

    :param source_path: The path to the file or folder to zip.
    :return: The path to the created zip file.
    """
    # Check if the source path is valid
    if not os.path.exists(source_path):
        raise FileNotFoundError(f"The path {source_path} does not exist.")

    if os.path.isdir(source_path):
        source_path = os.path.join(source_path, "")

    temp_dir = tempfile.mkdtemp()
    zip_path = os.path.join(temp_dir, constants.TMP_ZIP_FILENAME)

    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zipf:
        # If the source is a file, add it directly
        if os.path.isfile(source_path):
            zipf.write(source_path, os.path.basename(source_path))
        # If the source is a folder, walk through the folder and add files
        else:
            for root, _, files in os.walk(source_path):
                for file in files:
                    file_path = os.path.join(root, file)
                    zipf.write(
                        file_path,
                        os.path.relpath(file_path,
                                        os.path.dirname(source_path)))

    return zip_path


def _list_files(root_path: str, prefix_path: str) -> List[str]:
    file_list = []
    for dirpath, _, filenames in os.walk(root_path):
        for filename in filenames:
            relative_path = os.path.relpath(os.path.join(dirpath, filename),
                                            root_path)
            full_path = os.path.join(prefix_path, relative_path)
            file_list.append(full_path)
    return file_list

def remove_workspace(remote_dir, file_path=None) -> bool:
    """Removes a workspace folder or a workspace file.

    Args:
        remote_dir (str): The remote directory to remove.
        file_path (str, optional): The path of the file to remove. If not
            provided, the entire directory will be removed.
    
    Returns:
        True if the files were removed successfully, False otherwise.
    """
    api = storage_api.StorageApi(inductiva.api.get_client())

    logging.info("Removing workspace file(s)...")
    try:
        query_params = {"file_path": file_path} if file_path is not None else {}

        api.delete_file(
            query_params=query_params,
            path_params={
                "folder_name": remote_dir,
            },
        )
        logging.info("Workspace file(s) removed successfully.")
    except exceptions.ApiException as e:
        logging.error("An error occurred while removing the workspace:")
        logging.error(" > %s", json.loads(e.body)["detail"])
        return False
    return True
