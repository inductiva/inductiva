"""
Util functions to handle input and output data
when performing requests to the API.
"""
import os
import json
import zipfile
import tempfile
import shutil
import numpy as np

from absl import logging


def get_validate_request_params(original_params: dict,
                                type_annotations: dict) -> dict:
    """
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

    Return: dictionary with the params for request validation.
    """
    params = {}

    for variable in original_params:
        param_type = type_annotations[variable]

        if param_type == np.ndarray:
            params[variable] = {
                "shape": original_params[variable].shape,
            }
        else:
            params[variable] = original_params[variable]

    return params


def pack_param(name: str, value, param_type, dst_dir):
    """
    Pack a single parameter to the format used when passing the inputs
    to the web API.

    Args:
        name: Name of the parameter.
        value: Original value of the parameter.
        param_type: Type annotation of the param.
        dst_dir: Directory in which to store files that may be required for
            some param types.

    Return: Value that is passed to the API in a JSON file. For params
        that are sent as a file, the value is the path of the file relative
        to `dst_dir`.
    """
    if param_type == np.ndarray:
        param_filename = f"{name}.npy"
        param_fullpath = os.path.join(dst_dir, param_filename)
        np.save(param_fullpath, value)
        logging.debug("Stored %s to %s", name, param_fullpath)
        return param_filename

    return value


def pack_input(params, type_annotations, zip_name: str) -> str:
    """
    Pack all input params and compress all files into a zip file.
    All required files are created in a temporary directory which is then
    compressed with "zip" format. The path of the resulting zip file is
    returned.

    Args:
        params: Dict with the params that are passed into
            the request by the user.
        type_annotations: Dict with the type annotation of each param.
        zip_name: Name to use for the zip file (excluding the file extension).

    Return: Path to zip file with the compressed input. The zip will be located
        in the temporary directory of the OS (`/tmp` in linux).
    """
    with tempfile.TemporaryDirectory() as tmpdir_path:
        input_params = {}

        for variable in params:
            input_params[variable] = pack_param(
                name=variable,
                value=params[variable],
                param_type=type_annotations[variable],
                dst_dir=tmpdir_path,
            )

        # Write input dictionary with packed params to a JSON file
        input_json_path = os.path.join(tmpdir_path, "input.json")
        with open(input_json_path, "w", encoding="UTF-8") as fp:
            json.dump(input_params, fp)

        # Zip inputs
        zip_path = shutil.make_archive(
            os.path.join(tempfile.gettempdir(), zip_name), "zip", tmpdir_path)

        logging.debug("Compressed inputs to %s", zip_path)

    return zip_path


def unpack_value(value: str, var_type, output_dir: str):
    """
    Unpack a single output value, returning it as the correct type
    given the type annotation.

    Args:
        value: Output value to unpack. Can be a path relative to
            `output_dir` when the value is packed as a file.
        var_type: Type annotation of the value to unpack.
        output_dir: Directory where values packed as files are located.

    Return: Unpacked value with the type defined by `var_type`.
    """
    if var_type == np.ndarray:
        return np.load(os.path.join(output_dir, value))

    return type(value)


def unpack_output(zip_path: str, output_dir: str, return_type) -> any:
    """
    Unpack zip with the outputs of a task executed remotely in the API.

    Args:
        zip_path: Path to the zip file with the compressed outputs.
        output_dir: Directory in which to store the extracted outputs.
        return_type: Type annotation of the method return type. If `return_type`
            is a tuple, the return is also a tuple with the all output values
            unpacked.

    Return: Unpacked output of the method.
    """
    with zipfile.ZipFile(zip_path, "r") as zip_fp:
        zip_fp.extractall(output_dir)

    logging.debug("Extracted output to %s", output_dir)

    output_json_path = os.path.join(output_dir, "output.json")
    with open(output_json_path, "r", encoding="UTF-8") as fp:
        result_list = json.load(fp)

    if return_type.__name__ == "tuple":
        all_types = return_type.__args__

        return tuple(
            unpack_value(value, type, output_dir)
            for value, type in zip(result_list, all_types))

    return unpack_value(result_list[0], return_type, output_dir)
