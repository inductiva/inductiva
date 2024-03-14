"""Inductiva template manager class to manipulate files and render templates."""

from typing import Optional, Type
import pathlib
import shutil
import os

from inductiva import types
from inductiva.utils import file_manager
from inductiva.templating import renderers


class TemplateManager(file_manager.FileManager):
    """Template manager for rendering files and directories.

    This class extends functionality of the FileManager to render templates
    files. It provides methods to render single files and entire directories,
    """

    def __init__(self,
                 template_dir: types.PathOrStr,
                 root_dir: types.PathOrStr = "rendered_dir",
                 renderer_cls: Type[renderers.RendererAdapter] = None):
        """Initialize the template manager object.

        Args:
            template_dir (types.PathOrStr): Path to the template directory.
            root_dir (types.PathOrStr): Path to the default directory where
                rendered files will be written to.
            renderer (Type[renderers.RendererAdapter]): A Template Renderer
                class to use for rendering the templates. If not provided, the
                default renderer is used.
        """
        if not root_dir:
            raise ValueError("Root directory must be valid path or string")

        renderer_cls = renderer_cls or renderers.JinjaRendererAdapter
        self.set_root_dir(root_dir)
        self.renderer = renderer_cls(template_dir)
        self.template_dir = pathlib.Path(template_dir)
        self.template_dir_struct = self._gen_dir_struture()

    def _gen_dir_struture(self):
        """
        Generate a dictionary with the structure of the template directory.
        The keys are the subdirectories of the template directory, relative to
        the template directory itself. The values are dictionaries with the keys
        "files" and "templates" containing the lists of the files and templates
        in each subdirectory, respectively.
        Example:
        >>> self._gen_dir_struture()
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
        is_template = self.renderer.is_template
        return {
            os.path.relpath(root, start=self.template_dir): {
                "templates": [f for f in files if is_template(f)],
                "files": [f for f in files if not is_template(f)],
            } for root, _, files in os.walk(self.template_dir, topdown=True)
        }

    def render_file(self,
                    source_file: types.PathOrStr,
                    target_file: Optional[types.PathOrStr] = None,
                    overwrite: bool = False,
                    **render_args) -> pathlib.Path:
        """Render a single template to a file.

        Renders the `source_file` from the template directory to `target_file`
        inside the root directory, using the given render arguments.
        If `target_file` is not given, the name of `source_file` is used and the
        file is rendered directly to the root directory. When `overwrite` is
        False (default) and the destination files already exists, a
        FileExistsError will be raised.

        Args:
            source_file (str): The relative path of the template file inside
                the template directory.
            target_file: The name of the destination file. If not given, the
                name of the template will be used and the rendered file will be
                dumped to the root directory.
            overwrite: If True, the destination file will be overwritten if it
                already exists. If False (default), a FileExistsError will be
                raised if the destination file already exists.
            render_args: Keyword arguments that will be passed to the template
                renderer.

        Raises:
            FileExistsError: If the destination file already exists and
                `overwrite` is False.
        """

        if not self.renderer.is_template(source_file):
            raise ValueError(f"{source_file} is not a template file.")

        if target_file is None:
            target_file = pathlib.Path(source_file).name
        target_file = pathlib.Path(target_file)

        target_file = self.renderer.strip_extension(target_file)
        target_file = self.get_root_dir() / target_file

        if not overwrite and target_file.exists():
            raise FileExistsError(f"Target file {target_file} already exists.")

        target_file.parent.mkdir(parents=True, exist_ok=True)
        self.renderer.render(template_name=source_file,
                             target_file=target_file,
                             **render_args)

        return target_file

    def render_dir(self,
                   source_dir: Optional[types.PathOrStr] = None,
                   target_dir: Optional[types.PathOrStr] = None,
                   overwrite: bool = False,
                   **render_args) -> pathlib.Path:
        """
        Render an entire template directory.

        Renders the entire template directory to the given `target_dir`,
        relative to the root directory, using the provided rendering args.
        If `source_dir` is not provided (default), the entire template directory
        is rendered. If `target_dir` is provided, is is assumed it is relative
        to the template directory, in which case only that sub-directory is
        rendered. I any case, the tree structure of the `source_dir` is
        preserved in the `target_dir`.

        When `target_dir` does not exist, it will be created. If any destination
        file inside `target_dir` already exists, a FileExistsError will be
        raised when `overwrite` is False (default). Otherwise, files are
        overwritten.

        Args:
            source_dir (types.PathOrStr): Path to the source directory, relative
                to the template directory. If None, the full template directory
                is rendered.
            target_dir (types.PathOrStr): Path to the target directory, relative
                to the root directory. If None, the root directory is used.
                If set to an absolute path, the root directory is not used.
            overwrite (bool): If True, the destination files will be overwritten
                if they already exist. If False (default), a FileExistsError
                will be raised if any destination file already exists.
            render_args: Keyword arguments to render the template files.

        Raises:
            FileExistsError: If the destination file already exists and
                `overwrite` is False.
        """

        if source_dir is None:
            source_dir = "."
        source_dir = pathlib.Path(source_dir)

        if target_dir is None:
            target_dir = "."
        target_dir = pathlib.Path(target_dir)

        target_dir = self.get_root_dir() / target_dir

        if not overwrite:
            try:
                self._check_prerender_dir(target_dir)
            except FileExistsError as e:
                msg = f"{e}; set `overwrite=True` to overwrite existing files."
                raise FileExistsError(msg) from e

        for subdir, contents in self.template_dir_struct.items():
            target_subdir = target_dir / subdir
            target_subdir.mkdir(parents=True, exist_ok=True)
            source_dir = self.template_dir / subdir

            # simply copy non-template files to destination
            for file in contents["files"]:
                shutil.copy(source_dir / file, target_subdir)

            # render template files
            for file in contents["templates"]:
                target = target_subdir / self.renderer.strip_extension(file)
                template_name = (pathlib.Path(subdir) / file).as_posix()
                self.renderer.render(template_name, target, **render_args)

        return target_dir

    def _check_prerender_dir(self, target_dir: pathlib.Path):
        """Check if the destination filenames exist.

        Check if the destination filenames exist and raise a FileExistsError if
        any of them already exists. For template files, the extension is
        stripped before checking.

        Args:
            target_dir (pathlib.Path): The destination directory.
        """
        for subdir, contents in self.template_dir_struct.items():
            dest_subdir = target_dir / subdir
            for file in contents["files"] + contents["templates"]:
                target_name = self.renderer.strip_extension(dest_subdir / file)
                if target_name.exists():
                    raise FileExistsError(f"File {target_name} already exists.")
