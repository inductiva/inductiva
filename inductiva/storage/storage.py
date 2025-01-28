"""Methods to interact with the user storage resources."""

import logging
import math
import os
import pathlib
import time
import urllib
from enum import Enum
from typing import List, Literal, Optional, Tuple

import tqdm

import inductiva
from inductiva import constants
from inductiva.api import methods
from inductiva.client import exceptions, models
from inductiva.client.apis.tags import storage_api
from inductiva.utils import format_utils

_boto3_imported = True
try:
    import boto3

    logging.getLogger("botocore").setLevel(logging.WARNING)
except ImportError:
    _boto3_imported = False


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


class ExportDestination(Enum):
    AWS_S3 = "aws-s3"

    def __str__(self):
        return self.value


def initiate_multipart_upload(filename, bucket_name):
    """
    Initiate a multipart upload on S3 and return the UploadId.
    """
    s3_client = boto3.client("s3")
    response = s3_client.create_multipart_upload(
        Bucket=bucket_name,
        Key=filename,
    )
    return response["UploadId"]


def generate_presigned_url(upload_id, part_number, filename, bucket_name):
    """
    Generate a presigned URL for uploading a part to S3.
    """
    s3_client = boto3.client("s3")
    method_parameters = {
        "Bucket": bucket_name,
        "Key": filename,
        "PartNumber": part_number,
        "UploadId": upload_id,
    }

    return s3_client.generate_presigned_url(
        "upload_part",
        Params=method_parameters,
        ExpiresIn=3600,
    )


def generate_complete_multipart_upload_signed_url(
    upload_id,
    filename,
    bucket_name,
):
    """
    Generate a presigned URL for completing the multipart upload.
    """
    s3_client = boto3.client("s3")

    signed_url = s3_client.generate_presigned_url(
        ClientMethod="complete_multipart_upload",
        Params={
            "Bucket": bucket_name,
            "Key": filename,
            "UploadId": upload_id
        },
        HttpMethod="POST",
    )

    return signed_url


def get_file_size(file_path):
    api = storage_api.StorageApi(inductiva.api.get_client())

    contents = api.list_storage_contents({
        "path": file_path,
        "max_results": 2,
    }).body
    if len(contents) > 1:
        raise ValueError(f"Multiple files found at {file_path}. "
                         "Please specify a single file.")

    return list(contents.values())[0]["size_bytes"]


def get_multipart_parts(size, min_part_size=50 * 1024 * 1024):
    """
    Calculate the size of each part and the total number of parts

    The goal is to divide the data into parts `min_part_size` each:
    1. No more than 10,000 parts are created (maximum parts allowed by S3).
    2. The part size might be increased to avoid exceeding the part limit.

    Args:
        size (int): The total size of the file to be uploaded, in bytes.
        min_part_size (int): The minimum size of each part, in bytes.

    Returns:
        tuple: (part_size, part_count)
            - part_size (int): The size of each part in bytes.
            - part_count (int): The total number of parts.
    """
    max_parts = 10000

    if size <= min_part_size:
        return size, 1

    # Calculate the part size based on the smaller of two values:
    # - At least `min_part_size`
    # - Maximum size to ensure no more than 10,000 parts (size // 10000)
    max_allowed_part_size = size // max_parts
    part_size = max(min_part_size, max_allowed_part_size)

    part_count = math.ceil(size / part_size)

    return part_size, part_count


def multipart_upload(
    path,
    parts_size,
    upload_parts,
    complete_multipart_url,
):
    """
    Perform the multipart upload using the server.
    """
    api = storage_api.StorageApi(inductiva.api.get_client())

    api.export_multipart_files(
        body={
            "path": path,
            "parts_size": parts_size,
            "upload_parts": upload_parts,
            "complete_multipart_url": complete_multipart_url,
        })


def export_to_aws_s3(path_to_export, min_part_size_mb, filename, bucket_name):
    if not _boto3_imported:
        raise ValueError("boto3 is not installed. Please run "
                         "'pip install inductiva[aws]' to install it.")
    try:
        boto3.client("sts").get_caller_identity()
    except Exception as e:  # pylint: disable=broad-exception-caught
        raise (ValueError("AWS credentials not found. Please set your "
                          "AWS credentials with 'aws configure'.")) from e
    try:
        boto3.client("s3").head_bucket(Bucket=bucket_name)
    except Exception as e:  # pylint: disable=broad-exception-caught
        raise ValueError(f"Bucket {bucket_name} not found.") from e

    # Step 1: Get the file size
    file_size = get_file_size(path_to_export)

    # Step 2: Calculate the part size and count
    parts_size, parts_count = get_multipart_parts(
        file_size,
        min_part_size=min_part_size_mb * 1024 * 1024,
    )

    # Step 3: Initiate the multipart upload on aws
    upload_id = initiate_multipart_upload(filename, bucket_name)

    # Step 4: Generate presigned URLs for each part
    upload_parts = []
    for part_number in range(1, parts_count + 1):
        presigned_url = generate_presigned_url(upload_id, part_number, filename,
                                               bucket_name)
        upload_parts.append({
            "part_number": part_number,
            "part_url": presigned_url
        })

    # Step 5: Generate the complete multipart upload signed URL
    complete_multipart_url = generate_complete_multipart_upload_signed_url(
        upload_id,
        filename,
        bucket_name,
    )

    # Step 6: Ask the server to perform the multipart upload
    multipart_upload(
        path_to_export,
        parts_size,
        upload_parts,
        complete_multipart_url,
    )
    print(
        "Export is being done by inductiva server. You can close the terminal.")


def export(
    path_to_export: str,
    export_to: ExportDestination,
    bucket_name: str,
    file_name_to_save: str = None,
    min_part_size_mb: int = 50,
):
    file_name_to_save = file_name_to_save or pathlib.Path(path_to_export).name
    if export_to == ExportDestination.AWS_S3:
        print(f"Exporting {path_to_export} to {bucket_name}...")
        export_to_aws_s3(
            path_to_export,
            min_part_size_mb,
            file_name_to_save,
            bucket_name,
        )
    else:
        raise ValueError(f"Unsupported export destination: {export_to}")
