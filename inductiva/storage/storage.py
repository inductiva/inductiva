"""Methods to interact with the user storage resources."""
import json
import logging
import os
import tempfile
import zipfile
from typing import Literal
from urllib.parse import unquote, urlparse

import inductiva
from inductiva import constants
from inductiva.api import methods
from inductiva.client import exceptions
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
    file_name: str = None,
    unzip: bool = False,
):
    """
    Upload a file from a URL to the user workspace.

    Args:
        url (str): The URL of the file to upload.
        remote_dir (str): The remote directory to upload the file to. 
        file_name (str, optional): The name to save the file as. If not
            provided, the name will be extracted from the URL.
        unzip (bool, optional): Whether to unzip the file after uploading.
            Default is False.
    """
    api_instance = storage_api.StorageApi(inductiva.api.get_client())

    if file_name is None:
        parsed_url = urlparse(url)
        file_name = unquote(parsed_url.path.split("/")[-1])
    contents = api_instance.upload_from_url(
        query_params={
            "url": url,
            "file_name": file_name,
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
):
    """
    Upload a local file or directory to the user workspace.

    Args:
        local_path (str): The path to the local file or directory to be
            uploaded.
        remote_dir (str, optional): The remote directory where the file will
            be uploaded. Defaults to "default".
    """
    api_instance = storage_api.StorageApi(inductiva.api.get_client())

    input_zip_path = _zip_file_or_folder(local_path)

    zip_file_size = os.path.getsize(input_zip_path)
    logging.info("Input archive size: %s",
                 format_utils.bytes_formatter(zip_file_size))

    logging.info("Uploading input...")
    if methods.upload_file(
            api_instance=api_instance,
            input_zip_path=input_zip_path,
            remote_dir=remote_dir,
            get_upload_url_method=api_instance.get_upload_url,
            notify_upload_method=api_instance.notify_upload_file,
    ):
        logging.info("Input file successfully uploaded.")
    else:
        logging.error("An error occurred while uploading the input file.")

    logging.info("")
    os.remove(input_zip_path)


def _zip_file_or_folder(source_path):
    """
    Zips a file or a folder and saves it to a temporary folder.

    :param source_path: The path to the file or folder to zip.
    :return: The path to the created zip file.
    """
    # Check if the source path is valid
    if not os.path.exists(source_path):
        raise FileNotFoundError(f"The path {source_path} does not exist.")

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


def remove_workspace(remote_dir, file_name=None) -> bool:
    """Removes a workspace folder or a workspace file.

    Args:
        remote_dir (str): The remote directory to remove.
        file_name (str, optional): The name of the file to remove. If not
            provided, the entire directory will be removed.
    
    Returns:
        True if the files were removed successfully, False otherwise.
    """
    api = storage_api.StorageApi(inductiva.api.get_client())

    logging.info("Removing workspace file(s)...")
    try:
        query_params = {"file_name": file_name} if file_name is not None else {}

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
