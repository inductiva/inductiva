import os
import json
import zipfile
import tempfile
import shutil
import numpy as np

from absl import logging


def get_validate_request_params(original_params: dict, type_annotations):
    params = dict()

    for variable in original_params:
        param_type = type_annotations[variable]

        if param_type == np.ndarray:
            params[variable] = {
                "shape": original_params[variable].shape,
            }
        else:
            params[variable] = original_params[variable]

    return params


def pack_value(name, value, type, dir):
    if type == np.ndarray:
        param_filename = f"{name}.npy"
        param_fullpath = os.path.join(dir, param_filename)
        np.save(param_fullpath, value)
        logging.debug("Stored %s to %s", name, param_fullpath)
        return param_filename

    return value


def pack_input(params, type_annotations, zip_name: str) -> str:
    with tempfile.TemporaryDirectory() as tmpdir_path:
        input_params = dict()

        for variable in params:
            input_params[variable] = pack_value(
                name=variable,
                value=params[variable],
                type=type_annotations[variable],
                dir=tmpdir_path,
            )

        input_json_path = os.path.join(tmpdir_path, "input.json")
        with open(input_json_path, "w", encoding="UTF-8") as fp:
            json.dump(input_params, fp)

        zip_path = shutil.make_archive(
            os.path.join(tempfile.gettempdir(), zip_name), "zip", tmpdir_path)

        logging.debug("Compressed inputs to %s", zip_path)

    return zip_path


def unpack_value(value: str, type, output_path: str):
    if type == np.ndarray:
        return np.load(os.path.join(output_path, value))

    return type(value)


def unpack_output(zip_path: str, output_path: str, return_type) -> any:
    with zipfile.ZipFile(zip_path, "r") as zip_fp:
        zip_fp.extractall(output_path)

    logging.info("Extracted output to %s", output_path)

    with open(os.path.join(output_path, 'output.json'), "r") as fp:
        result_list = json.load(fp)

    if return_type.__name__ == "tuple":
        all_types = return_type.__args__

        return tuple(
            unpack_value(value, type, output_path)
            for value, type in zip(result_list, all_types))

    return unpack_value(result_list[0], return_type, output_path)
