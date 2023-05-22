"""Utils related to template files."""

from typing import Dict, List

from jinja2 import Environment, FileSystemLoader

from inductiva.utils.files import find_path_to_package

TEMPLATES_PATH = find_path_to_package("templates")


def replace_params_in_template(
    templates_dir: str,
    template_filename: str,
    params: Dict,
    output_file_path: str,
) -> None:
    """Replaces parameters in a template file."""

    environment = Environment(loader=FileSystemLoader(templates_dir))
    template = environment.get_template(template_filename)
    stream = template.stream(**params)
    stream.dump(output_file_path)


def batch_replace_params_in_template(
    templates_dir: str,
    template_filename_paths: List[str],
    params: Dict,
    output_filename_paths: List[str],
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

    for index, template_filename in enumerate(template_filename_paths):
        replace_params_in_template(templates_dir, template_filename, params,
                                   output_filename_paths[index])
