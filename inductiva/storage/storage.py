"""Methods to interact with the user storage resources."""
import datetime
import itertools
import logging
import math
import os
import pathlib
import threading
import time
import urllib
import zlib
from dataclasses import dataclass
from enum import Enum
from typing import List, Literal, Optional, Tuple

import concurrent.futures
import tqdm

import inductiva
from inductiva import constants
from inductiva import utils
from inductiva.api import methods
from inductiva.client import exceptions, models
from inductiva.client.apis.tags import storage_api
from inductiva.utils import format_utils
from inductiva.utils import data

MB = 1024 * 1024
_boto3_imported = True
try:
    import boto3
    from botocore.config import Config
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


def listdir(
    path="/",
    max_results: int = 10,
    order_by: Literal["size", "creation_time"] = "creation_time",
    sort_order: Literal["asc", "desc"] = "desc",
    print_results: bool = True,
):
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
    if len(path.split("/")) < 2:
        path += "/"

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
    if print_results:
        print(_print_contents_table(all_contents))

        _print_storage_size_and_cost()

        print(
            f"Listed {len(all_contents)} folder(s). Ordered by {order_by}.\n"
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


def get_signed_urls(
    paths: List[str],
    operation: Literal["upload", "download"],
) -> List[str]:
    api_instance = storage_api.StorageApi(inductiva.api.get_client())
    signed_urls = api_instance.get_signed_urls(query_params={
        "paths": paths,
        "operation": operation,
    }).body
    return signed_urls


@dataclass
class ZipFileInfo:
    """Represents information about a file within a ZIP archive."""
    name: str
    size: Optional[int]
    compressed_size: Optional[int]
    range_start: Optional[int]
    creation_time: Optional[datetime.datetime]
    compress_type: Optional[int]


@dataclass
class ZipArchiveInfo:
    """Represents the total ZIP size and file contents of a ZIP archive."""
    size: int
    files: List[ZipFileInfo]


def get_zip_contents(
    path: str,
    zip_relative_path: str = "",
    recursive: bool = False,
) -> ZipArchiveInfo:
    """
    Retrieve the contents of a ZIP archive from a given path.

    Args:
        path (str): The full path to the ZIP archive.
        zip_relative_path (str, optional): A relative path inside the ZIP 
            archive to filter the contents. Defaults to an empty string, 
            which lists all files within the archive.
        recursive (bool, optional): If True, list contents recursively within
            the specified `zip_relative_path`. If False, list only top-level 
            files and directories within the specified `zip_relative_path`.
            Defaults to False.

    Returns:
        ZipArchiveInfo: An object containing the total size of the ZIP archive
            and a list of `ZipFileInfo` objects representing the files 
            within the specified ZIP archive.
    """
    api_instance = storage_api.StorageApi(inductiva.api.get_client())
    query_params = {
        "path": path,
        "zip_relative_path": zip_relative_path,
        "recursive": str(recursive).lower()
    }
    response_body = api_instance.get_zip_contents(query_params).body
    files = [ZipFileInfo(
        name=str(file["name"]),
        size=int(file["size"]) \
            if file["size"] else None,
        compressed_size=int(file["compressed_size"]) \
            if file["compressed_size"] else None,
        range_start=int(file["range_start"]) \
            if file["range_start"] else None,
        creation_time=datetime.datetime.fromisoformat(file["creation_time"])
            if file["creation_time"] else None,
        compress_type=int(file["compress_type"])
            if file["compress_type"] else None
    ) for file in response_body["contents"]]
    return ZipArchiveInfo(size=int(response_body["size"]), files=files)


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
    },)
    logging.info("File is being uploaded...")
    logging.info("You can use 'inductiva storage ls' to check the status.")


def _convert_path_unix(path):
    """
    Convert a Windows path to a Unix-style path.
    """
    _, path = os.path.splitdrive(path)
    if "\\" in path:
        path = str(pathlib.PureWindowsPath(path).as_posix())
    return path


def _construct_remote_paths(local_path, remote_dir):
    """ Constructs remote paths for files to be uploaded."""
    remote_dir = os.path.normpath(remote_dir)
    local_path = os.path.normpath(_convert_path_unix(local_path))

    is_dir = os.path.isdir(local_path)
    if is_dir:
        local_dir = os.path.normpath(local_path)
        file_paths, total_size = _list_files(local_path)

        remote_file_paths = [
            os.path.normpath(os.path.join(remote_dir, file_path))
            for file_path in file_paths
        ]
    else:
        local_dir = os.path.dirname(local_path)
        filename = os.path.basename(local_path)
        remote_file_paths = [os.path.join(remote_dir, filename)]
        total_size = os.path.getsize(local_path)

    remote_file_paths = [
        _convert_path_unix(remote_path) for remote_path in remote_file_paths
    ]
    return remote_file_paths, local_dir, total_size


def upload(
    local_path: str,
    remote_dir: str,
):
    """
    Upload a local file or directory to a specified remote directory.

    Args:
        local_path (str): The path to the local file or directory to be
            uploaded.
        remote_dir (str): The remote directory where the file will
            be uploaded.
    
    Example:
        Upload a file to a remote directory:

        .. code-block:: python

            inductiva.storage.upload('local/path/file.txt', 'my_data')

        Upload a directory to a remote location:

        .. code-block:: python

            inductiva.storage.upload('local/path/folder', 'my_data')
    """

    if not os.path.exists(local_path):
        raise ValueError(f"File or directory '{local_path}' does not exist.")

    remote_file_paths, local_dir, total_size = _construct_remote_paths(
        local_path, remote_dir)

    if os.path.join(remote_dir, constants.TASK_OUTPUT_ZIP) in remote_file_paths:
        raise ValueError(f"Invalid file name: '{constants.TASK_OUTPUT_ZIP}.'")

    logging.info("Uploading content...")

    api_instance = storage_api.StorageApi(inductiva.api.get_client())

    urls = get_signed_urls(paths=remote_file_paths, operation="upload")

    with tqdm.tqdm(total=total_size,
                   unit="B",
                   unit_scale=True,
                   unit_divisor=1000) as progress_bar:

        for url, remote_file_path in zip(urls, remote_file_paths):
            file_path = remote_file_path.removeprefix(f"{remote_dir}/")
            local_file_path = _convert_path_unix(
                os.path.join(local_dir, file_path))

            try:
                methods.upload_file(api_instance, local_file_path, "PUT", url,
                                    progress_bar)
            except exceptions.ApiException as e:
                raise e

    logging.info("Content uploaded successfully.")


def _resolve_local_path(
    url,
    remote_base_path,
    local_base_dir,
    append_path=None,
    strip_zip=False,
):
    remote_url_path = urllib.parse.urlparse(url).path
    index = remote_url_path.find(remote_base_path)
    relative_path = remote_url_path[index:]

    if strip_zip and relative_path.endswith(".zip"):
        relative_path = relative_path.removesuffix(".zip")

    local_path = os.path.join(local_base_dir, relative_path)
    if append_path:
        local_path = os.path.join(local_path, append_path)

    os.makedirs(os.path.dirname(local_path), exist_ok=True)

    return local_path


def _get_size(url, pool_manager):
    response = pool_manager.urlopen("HEAD", url)
    size = int(response.getheader("Content-Length"))
    response.release_conn()
    return size


def _download_file(url,
                   download_path,
                   progress_bar,
                   progress_bar_lock,
                   pool_manager,
                   range_start=None,
                   range_end=None):
    headers = {}
    is_range = range_start is not None and range_end is not None

    if is_range:
        headers["Range"] = f"bytes={range_start}-{range_end}"
    response = pool_manager.urlopen("GET",
                                    url,
                                    headers=headers,
                                    preload_content=False)

    with open(download_path, "wb") as file:
        for chunk in response.stream():
            file.write(chunk)
            with progress_bar_lock:
                progress_bar.update(len(chunk))
    response.release_conn()

    return download_path


def _is_file_inside_zip(path):
    parts = path.split(os.sep)
    return any(part.endswith(".zip") for part in parts[:-1])


def _get_progress_bar(desc, total_bytes):
    return tqdm.tqdm(
        total=total_bytes,
        unit="B",
        unit_scale=True,
        unit_divisor=1000,
        desc=desc,
    )


def _download_file_from_inside_zip(remote_path, local_dir, pool_manager):
    before, after = remote_path.split(".zip" + os.sep, 1)
    path = before + ".zip"
    zip_relative_path = os.path.dirname(after)
    zip_filename = os.path.basename(after)

    url = get_signed_urls(paths=[path], operation="download")[0]

    if constants.TASK_OUTPUT_ZIP in path:
        prefix = data.ARTIFACTS_DIRNAME
    elif constants.TASK_INPUT_ZIP in path:
        prefix = data.INPUT_DIRNAME
    else:
        prefix = ""

    zip_files = get_zip_contents(path, prefix + zip_relative_path).files

    for zip_file in zip_files:
        if zip_file.name == zip_filename:
            break
    else:
        raise ValueError(f"File \"{after}\" not found in \"{path}\".")

    range_start = zip_file.range_start
    range_end = range_start + zip_file.compressed_size - 1
    compress_type = zip_file.compress_type
    file_path = after + ".zip" if compress_type else after
    download_path = _resolve_local_path(url, path, local_dir, file_path, True)

    desc = f"Downloading \"{zip_filename}\" from \"{remote_path}\""
    with _get_progress_bar(desc, range_start - range_end) as progress_bar:
        progress_bar_lock = threading.Lock()
        return _download_file(
            url=url,
            download_path=download_path,
            progress_bar=progress_bar,
            progress_bar_lock=progress_bar_lock,
            pool_manager=pool_manager,
            range_start=range_start,
            range_end=range_end,
        )


def _download_path(remote_path, local_dir, pool_manager):
    urls = get_signed_urls(paths=[remote_path], operation="download")
    with concurrent.futures.ThreadPoolExecutor() as executor:
        total_bytes = sum(
            executor.map(_get_size, urls, itertools.repeat(pool_manager)))

    resolved_paths = [
        _resolve_local_path(url, remote_path, local_dir) for url in urls
    ]

    desc = f"Downloading {len(urls)} file(s) from \"{remote_path}\""
    with _get_progress_bar(desc, total_bytes) as progress_bar:
        with concurrent.futures.ThreadPoolExecutor() as executor:
            progress_bar_lock = threading.Lock()
            return list(
                executor.map(_download_file, urls, resolved_paths,
                             itertools.repeat(progress_bar),
                             itertools.repeat(progress_bar_lock),
                             itertools.repeat(pool_manager)))


def _decompress(paths):
    for path in paths:
        decompress_dir, ext = os.path.splitext(path)
        # TODO: Improve the check for ZIP file
        if ext != ".zip":
            return
        utils.data.decompress_zip(path, decompress_dir)
        os.remove(path)


def _decompress_file_inside_zip(path):
    if not path.endswith(".zip"):
        return path

    with open(path, "rb") as f:
        compressed_data = f.read()
    decompressed_data = zlib.decompress(compressed_data, -zlib.MAX_WBITS)

    decompress_path = path.removesuffix(".zip")
    with open(decompress_path, "wb") as f:
        f.write(decompressed_data)
    os.remove(path)

    return decompress_path


def download(remote_path: str, local_dir: str = "", decompress: bool = True):
    """
    Downloads a file or folder from storage to a local directory, optionally 
    decompressing the contents.

    Args:
        remote_path (str): The path of the file or folder on the remote server
            to download.
        local_dir (str, optional): The local directory where the file or folder
            will be saved. Defaults to the current working directory.
        decompress (bool, optional): Whether to decompress the downloaded file 
            or folder if it is compressed. Defaults to True.

    Examples:
        Download a folder from a remote server to the current directory:

        .. code-block:: python

            inductiva.storage.download(remote_path="/path/to/remote/folder/")

        Download a file and save it to a local directory without decompressing:

        .. code-block:: python

            inductiva.storage.download(remote_path="/path/to/remote/file.zip",
                                       local_dir="/local/directory",
                                       decompress=False)

        Download a file inside a zip archive:

        .. code-block:: python

            inductiva.storage.download(
                remote_path="/some_task_id/output.zip/stdout.txt"
            )
    Note:
        It is not possible to download folders that are inside zip archives.
    """

    api_instance = storage_api.StorageApi(inductiva.api.get_client())
    pool_manager = api_instance.api_client.rest_client.pool_manager

    if _is_file_inside_zip(remote_path):
        download_path = _download_file_from_inside_zip(remote_path, local_dir,
                                                       pool_manager)
        paths = [_decompress_file_inside_zip(download_path)]
    else:
        paths = _download_path(remote_path, local_dir, pool_manager)

    if decompress:
        _decompress(paths)

    download_path = (paths[0] if len(paths) == 1 else os.path.join(
        local_dir, remote_path))

    logging.info('Successfully downloaded %d file(s) to "%s".', len(paths),
                 download_path)


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


def copy(source: str, target: str):
    """
    Copies a file or folder from a source path in storage to a target path.

    Args:
        source (str): The source path of the file or directory to copy.
        target (str): The destination path where the file or directory 
                      should be copied to.
    """
    api = storage_api.StorageApi(inductiva.api.get_client())
    api.copy(query_params={"source": source, "target": target})
    logging.info("Copied %s to %s successfully.", source, target)


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


def _initiate_multipart_upload(filename, bucket_name, region_name):
    """
    Initiate a multipart upload on S3 and return the UploadId.
    """
    s3_client = boto3.client("s3",
                             region_name=region_name,
                             config=Config(region_name=region_name,
                                           signature_version="v4"))
    response = s3_client.create_multipart_upload(
        Bucket=bucket_name,
        Key=filename,
    )
    return response["UploadId"]


def _generate_presigned_url(
    upload_id,
    part_number,
    filename,
    bucket_name,
    region_name,
):
    """
    Generate a presigned URL for uploading a part to S3.
    """
    s3_client = boto3.client("s3",
                             region_name=region_name,
                             config=Config(region_name=region_name,
                                           signature_version="v4"))
    method_parameters = {
        "Bucket": bucket_name,
        "Key": filename,
        "PartNumber": part_number,
        "UploadId": upload_id,
    }

    signed_url = s3_client.generate_presigned_url(
        "upload_part",
        Params=method_parameters,
        ExpiresIn=3600,
    )

    return signed_url


def _generate_complete_multipart_upload_signed_url(
    upload_id,
    filename,
    bucket_name,
    region_name,
):
    """
    Generate a presigned URL for completing the multipart upload.
    """
    s3_client = boto3.client("s3",
                             region_name=region_name,
                             config=Config(region_name=region_name,
                                           signature_version="v4"))

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


def _get_file_size(file_path):
    api = storage_api.StorageApi(inductiva.api.get_client())

    contents = api.list_storage_contents({
        "path": file_path,
        "max_results": 2,
    }).body
    if len(contents) > 1:
        raise ValueError(f"Multiple files found at {file_path}. "
                         "Please specify a single file.")

    return list(contents.values())[0]["size_bytes"]


def _get_multipart_parts(size: int,
                         part_size: int = 128 * MB) -> Tuple[int, int]:
    """
    Calculate the size of each part and the total number of parts

    The goal is to divide the data into parts `part_size` each:
    1. No more than 10,000 parts are created (maximum parts allowed by S3).
    2. The part size might be increased to avoid exceeding the part limit.
    3. The part size cannot be lower than 5MB.


    Args:
        size: The total size of the file to be uploaded, in bytes.
        part_size: The minimum size of each part, in bytes.

    Returns:
        - part_size: The size of each part in bytes.
        - part_count: The total number of parts.
    """
    max_parts = 10000
    min_allowed_part_size = 5 * MB  # 5MB

    # Ensure part_size is at least 5MB
    part_size = max(part_size, min_allowed_part_size)

    if size <= part_size:
        return size, 1

    # Calculate the part size based on the smaller of two values:
    # - At least `part_size`
    # - Maximum size to ensure no more than 10,000 parts (size // 10000)
    max_allowed_part_size = size // max_parts
    part_size = max(part_size, max_allowed_part_size)

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


def export_to_aws_s3(path_to_export, part_size, filename, bucket_name):
    if not _boto3_imported:
        print("boto3 is not installed. Please run "
              "'pip install inductiva[aws]' to install it.")
        return
    try:
        boto3.client("sts").get_caller_identity()
    except Exception:  # pylint: disable=broad-exception-caught
        print("AWS credentials not found. Please set your "
              "AWS credentials with 'aws configure'.")
        return
    try:
        boto3.client("s3").head_bucket(Bucket=bucket_name)
    except Exception:  # pylint: disable=broad-exception-caught
        print(f"Bucket {bucket_name} not found. Make sure the bucket exists "
              "and you have the correct permissions: "
              "https://tutorials.inductiva.ai/how_to/export-files-aws.html")
        return

    region_name = boto3.Session().region_name
    if not region_name:
        print("AWS region not found. Please set your AWS region with "
              "'aws configure'.")
        return

    # Step 1: Get the file size
    file_size = _get_file_size(path_to_export)

    # Step 2: Calculate the part size and count
    parts_size, parts_count = _get_multipart_parts(
        file_size,
        part_size=part_size * MB,
    )

    # Step 3: Initiate the multipart upload on aws
    upload_id = _initiate_multipart_upload(filename, bucket_name, region_name)

    # Step 4: Generate presigned URLs for each part
    upload_parts = []
    for part_number in range(1, parts_count + 1):
        presigned_url = _generate_presigned_url(upload_id, part_number,
                                                filename, bucket_name,
                                                region_name)
        upload_parts.append({
            "part_number": part_number,
            "part_url": presigned_url
        })

    # Step 5: Generate the complete multipart upload signed URL
    complete_multipart_url = _generate_complete_multipart_upload_signed_url(
        upload_id, filename, bucket_name, region_name)

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
    file_name: Optional[str] = None,
    part_size: int = 128,
):
    file_name = file_name or pathlib.Path(path_to_export).name
    if export_to == ExportDestination.AWS_S3:
        print(f"Exporting {path_to_export} to {bucket_name}...")
        export_to_aws_s3(
            path_to_export,
            part_size,
            file_name,
            bucket_name,
        )
    else:
        raise ValueError(f"Unsupported export destination: {export_to}")
