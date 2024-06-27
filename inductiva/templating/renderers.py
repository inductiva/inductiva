"""Jinja2-based renderer to manage template files."""
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

    def render_file(self, template_name: types.PathOrStr,
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

    def get_dir_structure(self, template_dir: types.PathOrStr):
        """
        Generate a dictionary with the structure of the template directory.
        The keys are the subdirectories of the template directory, relative to
        the template directory itself. The values are dictionaries with the keys
        "files" and "templates" containing the lists of the files and templates
        in each subdirectory, respectively.
        Example:
        >>> renderer.get_dir_structure("template_dir")
        {
            ".": {
                "files": ["file0.1"],
                "templates": ["template0.1"]
            },
            "subdir1": {
                "files": ["file1.1", "file1.2"],
                "templates": ["template1.1", "template1.2"]
            }
            ...
        }
        """
        return {
            os.path.relpath(root, start=template_dir): {
                "templates": [f for f in files if self.is_template(f)],
                "files": [f for f in files if not self.is_template(f)],
            } for root, _, files in os.walk(template_dir, topdown=True)
        }

    def check_prerender_dir(self, source_dir_struct, target_dir: pathlib.Path):
        """Check if the destination filenames exist.

        Check if the destination filenames exist and raise a FileExistsError if
        any of them already exists. For template files, the extension is
        stripped before checking.

        Args:
            source_dir_struct (dict): The structure of the source directory.
            target_dir (pathlib.Path): The destination directory.
        """
        for subdir, contents in source_dir_struct.items():
            dest_subdir = target_dir / subdir
            for file in contents["files"] + contents["templates"]:
                target_name = self.strip_extension(dest_subdir / file)
                if target_name.exists():
                    raise FileExistsError(f"File {target_name} already exists.")
