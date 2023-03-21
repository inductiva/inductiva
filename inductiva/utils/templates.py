"""Utils related to template files.

NB: These are temporary, and should be removed when we transition to a template
engine like jinja2.
"""

from typing import Dict

from jinja2 import Environment, FileSystemLoader


def replace_params_in_template_file(
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
