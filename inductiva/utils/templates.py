"""Utils related to template files."""

import os
import io
from typing import Dict, List, Union
from pathlib import Path

from jinja2 import Environment, FileSystemLoader

from inductiva.utils.files import find_path_to_package

TEMPLATES_PATH = find_path_to_package("templates")


def replace_params(
    template_path: str,
    params: Dict,
    output_file: Union[str, io.StringIO],
    remove_template: bool = False,
) -> None:
    """Replaces parameters in a template file."""

    template_path = Path(template_path)

    template_dir = template_path.parent
    template_filename = template_path.name

    environment = Environment(loader=FileSystemLoader(template_dir))
    template = environment.get_template(template_filename)
    stream = template.stream(**params)
    stream.dump(output_file)

    if remove_template:
        os.remove(template_path)


def batch_replace_params(
    templates_dir: str,
    template_filenames: List[str],
    params: Dict,
    output_filename_paths: List[str],
    remove_templates: bool = False,
) -> None:
    """Replaces parameters in a set of template files.

    For some simulators, more than one file needs to be changed.
    Moreover, some parameters are altered in different input files.
    To simplify, we can do a batch change for the same params dict.

    Args:
        templates_dir: Directory with the template files.
        template_filename_paths: List containing all template files to
            be changed.
        params: Dictionary of params that are inserted in the template
            files.
        output_filename_paths: List containing the output files in regard of
            the template files.
    """

    for index, template_filename in enumerate(template_filenames):
        template_path = os.path.join(templates_dir, template_filename)

        replace_params(
            template_path=template_path,
            params=params,
            output_file=output_filename_paths[index],
            remove_template=remove_templates,
        )
