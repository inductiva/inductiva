"""Inductiva template manager class to manipulate files and render templates."""
import pathlib
import shutil
from . import helpers

import jinja2


class TemplateManager:
    """Template manager for rendering files and directories.

    Template files contains variables that are replaced by
    values provided as keyword arguments to the render methods.
    The underlying template engine is the Jinja2 library.
    """

    @classmethod
    def render_dir(cls,
                   source_dir: str,
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
            overwrite (bool): If True, the destination folder will first
                be deleted if it already exists.
                If False (default), a FileExistsError will be raised if
                any file already exists in the target directory.
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

        loader = jinja2.FileSystemLoader(source_dir)
        environment = jinja2.Environment(loader=loader,
                                         keep_trailing_newline=True,
                                         undefined=jinja2.StrictUndefined)

        template_dir_struct = helpers.get_dir_structure(source_dir)

        if target_dir.exists() and any(target_dir.iterdir()):
            if overwrite:
                shutil.rmtree(target_dir)
            else:
                msg = (f"Non-empty target directory {target_dir}.\n"
                       "Use `overwrite=True` to delete previous content.")
                raise FileExistsError(msg)

        for subdir, contents in template_dir_struct.items():
            target_subdir = target_dir / subdir
            target_subdir.mkdir(parents=True, exist_ok=True)
            source_subdir = source_dir / subdir

            # simply copy non-template files to destination
            for file in contents["files"]:
                shutil.copy(source_subdir / file, target_subdir)

            # render template files
            for file in contents["templates"]:
                target_path = target_subdir / helpers.strip_extension(file)
                template_name = (pathlib.Path(subdir) / file).as_posix()
                template = environment.get_template(template_name)
                stream = template.stream(**render_args)
                stream.dump(str(target_path))
