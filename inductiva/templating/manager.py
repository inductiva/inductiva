import pathlib
import shutil
import os
from inductiva import types
from inductiva.templating import renderers
"""Inductiva template manager class to manipulate files and render templates."""


class TemplateManager:
    """Template manager for rendering files and directories.

    Template files contains variables that are replaced by
    values provided as keyword arguments to the render methods.
    The underlying template engine is the Jinja2 library.
    """

    @classmethod
    def render_dir(cls,
                   source_dir: types.PathOrStr,
                   target_dir: types.PathOrStr,
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
            source_dir (types.PathOrStr): Path to the source directory,
                containing template files (and potentially static files to be
                copied without modification).
            target_dir (types.PathOrStr): Path to the target directory.
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
        if not target_dir:
            raise ValueError("Target directory not provided.")

        source_dir = pathlib.Path(source_dir)
        target_dir = pathlib.Path(target_dir)

        renderer = renderers.JinjaRenderer(source_dir)

        template_dir_struct = renderer.get_dir_structure(source_dir)

        if not overwrite:
            try:
                renderer.check_prerender_dir(template_dir_struct, target_dir)
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
                target = target_subdir / renderer.strip_extension(file)
                template_name = (pathlib.Path(subdir) / file).as_posix()
                renderer.render_file(template_name, target, **render_args)

        return target_dir
