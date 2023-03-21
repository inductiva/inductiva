"""Utils related to template files.

NB: These are temporary, and should be removed when we transition to a template
engine like jinja2.
"""

import re
from typing import Dict


def replace_params_in_template_file(
    template_file_path: str,
    params: Dict,
    output_file_path: str,
) -> None:
    """Replaces parameters in a template file."""

    with open(template_file_path, "r", encoding="utf-8") as template_file:
        template_str = template_file.read()

    for key, value in params.items():
        template_str = re.sub(key, str(value), template_str)

    with open(output_file_path, "w", encoding="utf-8") as output_file:
        output_file.write(template_str)
