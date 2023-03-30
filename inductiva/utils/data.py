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
import numpy as np
import scipy

from absl import logging

from .meta import is_tuple
from inductiva.types import Path

INPUT_FILENAME = "input.json"
OUTPUT_FILENAME = "output.json"
ARTIFACTS_DIRNAME = "artifacts"


def get_validate_request_params(original_params: dict,
                                type_annotations: dict) -> dict:
    """Convert original params to the params used for request validation.

    Convert a dictionary with the request params into the params that
    are used to validate the request.
    For instance, a numpy array is replaced by its shape,
    as that is what's necessary to validate if the user has permissions
    to execute the request.
    The original params are not modified.

    Args:
        original_params: Dict with the params that are passed into
            the request by the user.
        type_annotations: Dict with the type annotation of each param.

    Return:
        Returns a dictionary with the params for request validation.
    """
    params = {}

    for variable in original_params:
        param_type = type_annotations.get(variable, None)
        if param_type in (np.ndarray, scipy.sparse):
            params[variable] = {
                "shape": original_params[variable].shape,
            }
        elif param_type == pathlib.Path:
            # TODO: what kind of information do we want to send for validation
            # of a SplishSplash simulation?
            params[variable] = str(original_params[variable])
        else:
            params[variable] = original_params[variable]

    return params


def pack_param(name: str, value, param_type, dst_dir):
    """Pack a single parameter to the format expected by the Web API.

    Args:
        name: Name of the parameter.
        value: Original value of the parameter.
        param_type: Type annotation of the param.
        dst_dir: Directory in which to store files that may be required for
            some param types.

    Return:
        Returns a value that is passed to the API in a JSON file. For params
        that are sent as a file, the value is the path of the file relative
        to `dst_dir`.
    """
    if param_type == np.ndarray:
        param_filename = f"{name}.npy"
        param_fullpath = os.path.join(dst_dir, param_filename)
        np.save(param_fullpath, value)
        logging.debug("Stored %s to %s", name, param_fullpath)
        return param_filename

    if param_type == scipy.sparse:
        param_filename = f"{name}.npz"
        param_fullpath = os.path.join(dst_dir, param_filename)
        scipy.sparse.save_npz(param_fullpath, value)
        logging.debug("Stored %s to %s", name, param_fullpath)
        return param_filename

    if param_type == pathlib.Path:
        dst_dir_name = name
        dst_fullpath = os.path.join(dst_dir, dst_dir_name)

        shutil.copytree(value, dst_fullpath)

        logging.debug("Copied %s to %s", value, dst_fullpath)
        return str(dst_dir_name)

    return value


def pack_input(params, type_annotations, zip_name) -> str:
    """Pack all inputs into a zip file.

    Pack all input params and compress all files into a zip file.
    All required files are created in a temporary directory which is then
    compressed with "zip" format. The path of the resulting zip file is
    returned.

    Args:
        params: Dict with the params that are passed into
            the request by the user.
        type_annotations: Dict with the type annotation of each param.
            If a type annotation doesn't exist for a param, it is
            assumed that it is JSON encodable.
        zip_name: Name of the zip file to be created.

    Return:
        Returns a path to zip file with the compressed input. The zip will be
        located in the temporary directory of the OS (`/tmp` in linux).
    """
    with tempfile.TemporaryDirectory() as tmpdir_path:
        input_params = {}

        for variable in params:
            input_params[variable] = pack_param(
                name=variable,
                value=params[variable],
                param_type=type_annotations.get(variable, None),
                dst_dir=tmpdir_path,
            )

        # Write input dictionary with packed params to a JSON file
        input_json_path = os.path.join(tmpdir_path, INPUT_FILENAME)
        with open(input_json_path, "w", encoding="UTF-8") as fp:
            json.dump(input_params, fp)

        # Zip inputs
        zip_path = shutil.make_archive(
            os.path.join(tempfile.gettempdir(), zip_name), "zip", tmpdir_path)

        logging.debug("Compressed inputs to %s", zip_path)

    return zip_path


def unpack_value(value: str, var_type, output_dir: Path):
    """Unpack a single output value return by the Web API.

    Unpack a single output value, returning it as the correct type
    given the type annotation.

    Args:
        value: Output value to unpack. Can be a path relative to
            `output_dir` when the value is packed as a file.
        var_type: Type annotation of the value to unpack.
        output_dir: Directory where values packed as files are located.

    Return:
        Returns the unpacked value with the type defined by `var_type`.
    """
    if var_type == np.ndarray:
        return np.load(os.path.join(output_dir, value))

    if var_type == pathlib.Path:
        return pathlib.Path(os.path.join(output_dir, value))

    return type(value)


def unpack_output(result_list, output_dir: Path, return_type):
    """Transform outputs of a task executed in the API to the right objects.

    Args:
        result_list: List
        output_dir: Directory in which to store the extracted outputs.
        return_type: Type annotation of the method return type. If `return_type`
            is a tuple, the return is also a tuple with the all output values
            unpacked.

    Return:
        Returns the unpacked output of the method.
    """
    if return_type is None:
        return

    if return_type == pathlib.Path:
        return pathlib.Path(output_dir)

    if is_tuple(return_type):
        all_types = return_type.__args__

        return tuple(
            unpack_value(value, type, output_dir)
            for value, type in zip(result_list, all_types))

    return unpack_value(result_list[0], return_type, output_dir)


def extract_output(zip_path, output_dir):
    with zipfile.ZipFile(zip_path, "r") as zip_fp:
        output_json = zip_fp.read(OUTPUT_FILENAME)
        result_list = json.loads(output_json)

        extract_subdir_files(zip_fp, ARTIFACTS_DIRNAME, output_dir)
    return result_list


def extract_subdir_files(zip_fp: zipfile.ZipFile, dir_name: str,
                         output_dir: Path):
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
