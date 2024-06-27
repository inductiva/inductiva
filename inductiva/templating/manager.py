"""Inductiva template manager class to manipulate files and render templates."""
from typing import Union, IO
import pathlib
import shutil
import os

import jinja2


class TemplateManager:
    """Template manager for rendering files and directories.

    Template files contains variables that are replaced by
    values provided as keyword arguments to the render methods.
    The underlying template engine is the Jinja2 library.

    Attributes:
        TEMPLATE_EXTENSION (str): The file extension used by Jinja2 templates.
    """

    TEMPLATE_EXTENSION = ".jinja"

    def __init__(self, template_dir: str):
        """Initialize the template manager.

        Args:
            template_dir (str): The directory where the templates are
                located.
        """
        self.template_dir = template_dir

        loader = jinja2.FileSystemLoader(self.template_dir)
        self.environment = jinja2.Environment(loader=loader,
                                              undefined=jinja2.StrictUndefined)

    def render_file(self, template_name: str, target_file: Union[str, IO],
                    **render_args):
        """Render a template to a file.

        The template with the given `template_name` is rendered to a file named
        `dest_name`, using the given render arguments. If the template is not
        found, a jinja2.exceptions.TemplateNotFound exception will be raised.
        If the file pointed by `template_name` is not a template, the file will
        be copied to `dest_name` as is.

        Args:
            template_name (str): The name of the template file to be rendered,
                relative to the template directory.
            dest_name (str): The name of the destination file.
            render_args: Keyword arguments that will be passed to the template
                renderer.
        """
        # jinja uses posix paths (forward slashes) to find templates,
        # even on Windows
        template_name = pathlib.Path(template_name).as_posix()
        template = self.environment.get_template(template_name)
        stream = template.stream(**render_args)
        stream.dump(target_file)

    @staticmethod
    def render_dir(source_dir: str,
                   target_dir: str,
                   overwrite: bool = False,
                   **render_args):
        """
        Renders the entire template directory to the given `target_dir`,
        using the provided rendering args, preserving the directory structure of
        the `source_dir`. Files that are not templates are copied as is.

        When `target_dir` does not exist, it will be created. If any destination
        file inside `target_dir` already exists, a FileExistsError will be
        raised when `overwrite` is False (default). Otherwise, files are
        overwritten.

        Args:
            source_dir (str): Path to the source directory,
                containing template files (and potentially static files to be
                copied without modification).
            target_dir (str): Path to the target directory.
            overwrite (bool): If True, the destination files will be overwritten
                if they already exist. If False (default), a FileExistsError
                will be raised if any destination file already exists.
            render_args: Keyword arguments to render the template files.

        Raises:
            FileExistsError: If the destination file already exists and
                `overwrite` is False.
        """
        if not source_dir:
            raise ValueError("Source directory not provided.")

        source_dir = pathlib.Path(source_dir)
        if not source_dir.is_dir():
            raise ValueError(f"Source directory {source_dir} does not exist.")

        if not target_dir:
            raise ValueError("Target directory not provided.")

        target_dir = pathlib.Path(target_dir)

        manager = TemplateManager(source_dir)

        template_dir_struct = get_dir_structure(source_dir)

        if not overwrite:
            try:
                check_prerender_dir(template_dir_struct, target_dir)
            except FileExistsError as e:
                msg = f"{e}; set `overwrite=True` to overwrite existing files."
                raise FileExistsError(msg) from e

        for subdir, contents in template_dir_struct.items():
            target_subdir = target_dir / subdir
            target_subdir.mkdir(parents=True, exist_ok=True)
            source_subdir = source_dir / subdir

            # simply copy non-template files to destination
            for file in contents["files"]:
                shutil.copy(source_subdir / file, target_subdir)

            # render template files
            for file in contents["templates"]:
                target = target_subdir / strip_extension(file)
                template_name = (pathlib.Path(subdir) / file).as_posix()
                manager.render_file(template_name, str(target), **render_args)


def is_template(file: str) -> bool:
    """Check if the given file is of template type."""
    return file.endswith(TemplateManager.TEMPLATE_EXTENSION)


def strip_extension(file: str) -> str:
    """Strip the template extension if the given file is a template."""
    if not is_template(file):
        return file
    return os.path.splitext(file)[0]


def get_dir_structure(template_dir: str):
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
            "templates": [f for f in files if is_template(f)],
            "files": [f for f in files if not is_template(f)],
        } for root, _, files in os.walk(template_dir, topdown=True)
    }


def check_prerender_dir(source_dir_struct, target_dir: str):
    """Check if the destination filenames exist.

    Check if the destination filenames exist and raise a FileExistsError if
    any of them already exists. For template files, the extension is
    stripped before checking.

    Args:
        source_dir_struct (dict): The structure of the source directory,
            as returned by `get_dir_structure`.
        target_dir (pathlib.Path): The destination directory.
    """
    target_dir = pathlib.Path(target_dir)
    for subdir, contents in source_dir_struct.items():
        dest_subdir = target_dir / subdir
        for file in contents["files"] + contents["templates"]:
            raw_target_name = str(dest_subdir / file)
            target_name = strip_extension(raw_target_name)
            if os.path.exists(target_name):
                raise FileExistsError(f"File {target_name} already exists.")
