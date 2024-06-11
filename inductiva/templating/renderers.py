"""Base and specific renderers packages to manage template files."""
from typing import Union, IO
import pathlib
import os

import jinja2

from inductiva import types


class JinjaRenderer():
    """Renderer based on Jinja2 template engine.

    Attributes:
        TEMPLATE_EXTENSION (str): The file extension used by Jinja2 templates.
    """

    TEMPLATE_EXTENSION = ".jinja"

    def __init__(self, template_dir: types.PathOrStr):
        """Initialize the Jinja based Renderer.
        Args:
            template_dir (PathOrStr): The directory where the templates are
                located.
        """
        self.template_dir = pathlib.Path(template_dir)

        loader = jinja2.FileSystemLoader(self.template_dir)
        self.environment = jinja2.Environment(loader=loader,
                                              undefined=jinja2.StrictUndefined)

    def render(self, template_name: types.PathOrStr,
               target_file: Union[types.PathOrStr, IO], **render_args):
        """Render a template to a file.

        The template with the given `template_name` is rendered to a file named
        `dest_name`, using the given render arguments. If the template is not
        found, a jinja2.exceptions.TemplateNotFound exception will be raised.
        If the file pointed by `template_name` is not a template, the file will
        be copied to `dest_name` as is.

        Args:
            template_name (str): The name of the template file to be rendered.
            dest_name (str or Path): The name of the destination file.
            render_args: Keyword arguments that will be passed to the template
                renderer.
        """
        # jinja uses posix paths (forward slashes) to find templates,
        # even on Windows
        template_name = pathlib.Path(template_name).as_posix()
        template = self.environment.get_template(template_name)
        stream = template.stream(**render_args)

        if isinstance(target_file, pathlib.Path):
            stream.dump(str(target_file))
        else:
            stream.dump(target_file)

    def is_template(self, file: types.PathOrStr) -> bool:
        """Check if the given file is of template type."""
        if isinstance(file, pathlib.Path):
            return file.suffix == self.TEMPLATE_EXTENSION
        elif isinstance(file, str):
            return file.endswith(self.TEMPLATE_EXTENSION)
        return False

    def strip_extension(self, file: types.PathOrStr) -> str:
        """Strip the template extension if the given file is a template."""
        if not self.is_template(file):
            return file
        if isinstance(file, pathlib.Path):
            return file.with_suffix("")
        return os.path.splitext(file)[0]
