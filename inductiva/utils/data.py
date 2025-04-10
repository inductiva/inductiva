"""Util functions to handle input and output data for API requests.

This module contains several functions related to packing and unpacking
of inputs and ouputs for Web API requests.
Additionally, it contains some global constant variables defining
configurations related to paths where certain files are expected to be.
"""
import os
import json
import pathlib
import zipfile
import tempfile
import shutil
from typing import Callable, List
from tqdm import tqdm
import fsspec
import urllib3

import logging

INPUT_FILENAME = "input.json"
OUTPUT_FILENAME = "output.json"
ARTIFACTS_DIRNAME = "artifacts"
INPUT_DIRNAME = "sim_dir"


def pack_input(input_dir, kwargs, zip_name) -> str:
    """Pack all inputs into a zip file.

    Pack all input params and compress all files into a zip file.
    All required files are created in a temporary directory which is then
    compressed with "zip" format. The path of the resulting zip file is
    returned.

    Args:
        input_dir: Directory containing the input files to be uploaded.
        kwargs: Dict with the params that are passed into
            the request by the user.
        zip_name: Name of the zip file to be created.

    Return:
        Returns a path to zip file with the compressed input. The zip will be
        located in the temporary directory of the OS (`/tmp` in linux).
    """
    with tempfile.TemporaryDirectory() as tmpdir_path:
        # Copy input directory to a temporary directory
        dst_fullpath = os.path.join(tmpdir_path, INPUT_DIRNAME)

        if input_dir:
            shutil.copytree(input_dir, dst_fullpath)

        # Write input dictionary with packed params to a JSON file
        input_json_path = os.path.join(tmpdir_path, INPUT_FILENAME)
        with open(input_json_path, "w", encoding="UTF-8") as fp:
            json.dump(kwargs, fp)

        # Zip everything in the temporary directory into a single zip file
        zip_path = shutil.make_archive(
            os.path.join(tempfile.gettempdir(), zip_name), "zip", tmpdir_path)

        logging.debug("Compressed inputs to %s", zip_path)

    return zip_path


def extract_output(zip_path, output_dir):
    with zipfile.ZipFile(zip_path, "r") as zip_fp:
        output_json = zip_fp.read(OUTPUT_FILENAME)
        result_list = json.loads(output_json)

        extract_subdir_files(zip_fp, ARTIFACTS_DIRNAME, output_dir)
    return result_list


def extract_subdir_files(zip_fp: zipfile.ZipFile, dir_name: str,
                         output_dir: pathlib.Path):
    """Util function to extract the contents of a directory in a ZIP archive.

    For instance, if a ZIP archive contains a directory called `dir_name`,
    the contents of that directory are extracted directly to `output_dir`.

    Args:
        zip_fp: ZipFile from which to extract the directory.
        dir_name: Name of the directory inside the ZIP archive.
        output_dir: Destination directory of the contents of `data_dir`.
    """
    for member in zip_fp.namelist():
        is_dir = not os.path.basename(member)

        if not member.startswith(dir_name) or is_dir:
            continue

        src_file = zip_fp.open(member)
        target_relative_path = pathlib.Path(member).relative_to(dir_name)
        target_path = os.path.join(output_dir, target_relative_path)

        os.makedirs(os.path.dirname(target_path), exist_ok=True)

        with open(target_path, "wb") as f:
            shutil.copyfileobj(src_file, f)


def zip_dir(dir_path, zip_name):
    """Compress a directory into a zip file."""
    zip_path = shutil.make_archive(
        os.path.join(tempfile.gettempdir(), zip_name), "zip", dir_path)

    logging.debug("Compressed inputs to %s", zip_path)

    return zip_path


def _extract_zip_file_to_dir(
    dest_dir: pathlib.Path,
    remove_zip_file: zipfile.ZipFile,
    filename: str,
    zip_path: str,
):
    """Write a file from a ZIP archive to the output directory.

    Args:
        dest_dir: Directory where to store the extracted file.
        remove_zip_file: ZipFile object from which to extract the file.
        filename: Name of the file to extract.
        zip_path: Path of the file inside the ZIP archive.
    """
    with remove_zip_file.open(zip_path) as source:
        target_path = dest_dir / pathlib.Path(filename)
        target_path.parent.mkdir(parents=True, exist_ok=True)
        with open(target_path, "wb") as target:
            target.write(source.read())


def _download_partial_files(
    download_url: str,
    filenames: List[str],
    dest_dir: pathlib.Path,
    make_zip_path: Callable,
) -> None:
    """Download the partial files of a task.

    Args:
        download_url: URL from which to download the files.
        filenames: List of filenames to download.
        dest_dir: Path where to store the downloaded files.
        make_zip_path: Function to create a path of a file inside the ZIP 
        given a filename.

    Return:
        Returns True if the download was successful, False otherwise.
    """
    try:
        remote_filesystem = fsspec.filesystem("http")

        with remote_filesystem.open(download_url, "rb") as remote_file:
            with zipfile.ZipFile(remote_file) as remote_zip_file:
                for filename in filenames:
                    zip_path = make_zip_path(filename)
                    try:
                        _extract_zip_file_to_dir(dest_dir, remote_zip_file,
                                                 filename, zip_path)
                    except KeyError:
                        logging.warning(
                            "File %s not found in the output archive.",
                            zip_path)

    except Exception as e:  # pylint: disable=broad-except
        logging.debug("Error downloading partial outputs: %s", e)
        logging.error("Partial download failed.")


def download_partial_outputs(
    download_url: str,
    filenames: List[str],
    output_dir: pathlib.Path,
) -> None:
    """Download the partial outputs of a task.

    Args:
        download_url: URL from which to download the outputs.
        filenames: List of filenames to download.
        output_dir: Path where to store the downloaded files.

    Return:
        Returns True if the download was successful, False otherwise.
    """
    return _download_partial_files(
        download_url=download_url,
        filenames=filenames,
        dest_dir=output_dir,
        make_zip_path=lambda filename: "artifacts/" + filename)


def download_partial_inputs(
    download_url: str,
    filenames: List[str],
    input_dir: pathlib.Path,
) -> None:
    """Download the partial inputs of a task.

    Args:
        download_url: URL from which to download the inputs.
        filenames: List of filenames to download.
        output_dir: Path where to store the downloaded files.

    Return:
        Returns True if the download was successful, False otherwise.
    """
    return _download_partial_files(
        download_url=download_url,
        filenames=filenames,
        dest_dir=input_dir,
        make_zip_path=lambda filename: "sim_dir/" + filename \
            if filename != "input.json" else filename)


def download_file(
    response: urllib3.response.HTTPResponse,
    output_path: pathlib.Path,
    chunk_size: int = 1000,
) -> None:
    """Download a file from a urllib3 response object.

    Use a urllib3 response object to download a file, showing a progress bar.

    Args:
        response: urllib3 response object.
        output_path: Path where to store the downloaded file.
        chunk_size: Size of the chunks in which to download the file.
    """
    # if the response header does not contain a x-content-length header,
    # the progress bar will not be displayed correctly, but download
    # will still work.
    # "x-content-length" is a custom header that is set by the API instead
    # of the standard "content-length" header, because the API needs to use
    # "transfer-encoding: chunked" to stream the response.
    # If the download is provided by a file server, the "content-length"
    # header will be used.
    download_size = response.headers.get(
        "x-content-length") or response.headers.get("content-length", 0)

    with tqdm(
            total=int(download_size),
            unit="B",
            unit_scale=True,
            unit_divisor=1000,  # Use 1 KB = 1000 bytes
    ) as progress_bar:
        with open(output_path, "wb") as f:
            while chunk := response.read(chunk_size):
                f.write(chunk)
                progress_bar.update(len(chunk))

    response.release_conn()


def decompress_zip(zip_path: pathlib.Path, output_dir: pathlib.Path):
    """Uncompress a ZIP archive containing the outputs of a task.

    If the archive contains the directory called artifacts, it means that the
    download includes the full outputs of the task with the full directory
    structure of the outputs (output.json, artifacts/*). Only the contents
    inside artifacts are extracted to the output directory. If the archive
    does not contain the artifacts directory, it means that the download
    only includes a few files, which are directly in the root of the archive,
    without the `artifacts` directory.
    """
    with zipfile.ZipFile(zip_path, "r") as zip_f:
        full_output = True
        try:
            zip_f.getinfo("artifacts/")
        except KeyError:
            full_output = False

        if full_output:
            extract_subdir_files(
                zip_f,
                "artifacts",
                output_dir,
            )
        else:
            zip_f.extractall(output_dir)
