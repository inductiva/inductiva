"""
Template rendering utilities.

Provides utility classes for rendering template files and directories.
"""
from contextlib import contextmanager
from abc import ABC, abstractmethod
from typing import Dict, Any, Union, IO
import pathlib
import shutil
import os

import jinja2

from inductiva.types import PathOrStr, OptionalPathOrStr
# PathOrStr = OptionalPathOrStr = Any

class RendererAdapter(ABC):
    """Abstract base class for template renderers."""

    @abstractmethod
    def render(self, template_name: str, dest_name: PathOrStr, **render_args):
        """Render a template to a file."""
        pass

    @abstractmethod
    def is_template(self, file: PathOrStr) -> bool:
        """Check if a file is a template."""
        pass

    @abstractmethod
    def set_template_dir(self, template_dir: PathOrStr) -> None:
        """Set the directory where the templates are located."""
        pass

    @abstractmethod
    def strip_extension(self, file: PathOrStr) -> str:
        """Strip the extension of a template file."""
        pass


class JinjaRendererAdapter(RendererAdapter):
    """Adapter for the Jinja2 template engine.

    Attributes:
        TEMPLATE_EXTENSION (str): The file extension used by Jinja2 templates.
    """

    TEMPLATE_EXTENSION = ".jinja"

    def __init__(self, template_dir: PathOrStr):
        """Initialize the JinjaRendererAdapter.

        Args:
            template_dir (PathOrStr): The directory where the templates are
                located.
        """
        self.set_template_dir(template_dir)

    def set_template_dir(self, template_dir: PathOrStr) -> None:
        """Set the directory where the templates are located.

        Sets the directory where the templates are located and initializes the
        Jinja2 environment and loader.

        Args:
            template_dir (PathOrStr): The directory where the templates are
                located.
        """
        self.template_dir = pathlib.Path(template_dir)
        self.loader = jinja2.FileSystemLoader(self.template_dir)
        self.environment = jinja2.Environment(loader=self.loader,
                                              undefined=jinja2.StrictUndefined)

    def render(self, template_name: str, dest_name: Union[PathOrStr, IO],
               **render_args) -> None:
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

        if isinstance(dest_name, pathlib.Path):
            stream.dump(str(dest_name))
        else:
            stream.dump(dest_name)


    @staticmethod
    def is_template(file: PathOrStr) -> bool:
        """Check if the given file is a Jinja2 template."""
        if isinstance(file, pathlib.Path):
            return file.suffix == JinjaRendererAdapter.TEMPLATE_EXTENSION
        return file.endswith(JinjaRendererAdapter.TEMPLATE_EXTENSION)

    @staticmethod
    def strip_extension(file: PathOrStr) -> str:
        """Strip the template extension if the given file is a template."""
        if not JinjaRendererAdapter.is_template(file):
            return file
        if isinstance(file, pathlib.Path):
            return file.with_suffix("")
        return os.path.splitext(file)[0]


class Template:
    """A class to render templates to files."""

    def __init__(self,
                 template_dir: PathOrStr,
                 renderer_cls: RendererAdapter = JinjaRendererAdapter):
        """Initialize the Template object.

        Args:
            template_dir (PathOrStr): The directory where the templates are
                located.
            renderer_cls (RendererAdapter): The class of the template renderer
                to be used. Defaults to JinjaRendererAdapter.
        """
        self.template_dir = pathlib.Path(template_dir)
        self.renderer = renderer_cls(self.template_dir)
        self.dir_struct = self._gen_dir_struture()
        self._per_template_ctx = {}

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

    @contextmanager
    def with_context(self, context: Dict[str, Dict[str, Any]]):
        """
        Define a per-template file context for rendering templates.

        The context is a dictionary with the template file names as keys and
        dictionaries with the rendering arguments as values. This context will
        be used to render the templates inside the context manager block when
        performing calls to the `render_dir` and `render_file` methods.
        Note that the rendering arguments given to the `render_dir` and
        `render_file` methods will overwrite the ones defined in the context.
        """
        ctx = self._per_template_ctx
        self._per_template_ctx = context
        try:
            yield self
        finally:
            self._per_template_ctx = ctx

    def render_dir(self,
                   dest_dir: OptionalPathOrStr = None,
                   overwrite: bool = False,
                   **render_args):
        """
        Render the entire template director to another directory.

        Renders the entire template directory to the given directory `dest_dir`.
        The directory structure is preserved and the files are rendered
        using the provided rendering args. If `dest_dir` is not provided
        (default), the current working directory is used.
        When `dest_dir` does not exist, it will be created. If any destination
        file inside `dest_dir` already exists, a FileExistsError will be
        raised when `overwrite` is False (default). Otherwise, files are
        overwritten. No check is made for the existence of `dest_dir` to
        prevent errors from being raised when files are to be rendered to the
        current working directory.

        Args:
            dest_dir (str or Path): The directory where the rendered files will
                be placed. If not given (the default), the current working
                directory will be used.
            overwrite: If True, the destination files will be overwritten if
                they already exist. If False (default), a FileExistsError will
                be raised if any destination file already exists.
            render_args: Keyword arguments that will be passed to the template
                renderer.
        """
        if dest_dir is None:
            dest_dir = pathlib.Path.cwd()
        dest_dir = pathlib.Path(dest_dir)

        context = self._per_template_ctx
        if not overwrite:
            try:
                self._validate_destination(dest_dir)
            except FileExistsError as e:
                msg = f"{e}; set `overwrite=True` to overwrite existing files."
                raise FileExistsError(msg) from e

        for subdir, contents in self.dir_struct.items():
            dest_subdir = dest_dir / subdir
            dest_subdir.mkdir(parents=True, exist_ok=True)
            source_dir = self.template_dir / subdir

            # simply copy non-template files to destination
            for file in contents["files"]:
                shutil.copy(source_dir / file, dest_subdir)

            # render template files
            for file in contents["templates"]:
                dest_name = dest_subdir / self.renderer.strip_extension(file)
                template_name = (pathlib.Path(subdir) / file).as_posix()
                self.renderer.render(template_name, dest_name,
                                     **context.get(template_name, {}),
                                     **render_args)

    def render_file(self,
                    template_name: str,
                    dest_name: OptionalPathOrStr = None,
                    dest_dir: OptionalPathOrStr = None,
                    overwrite: bool = False,
                    **render_args):
        """
        Render a template to a file.

        Renders the template with the given `template_name` to a file named
        `dest_name` inside the directory `dest_dir`, using the given render
        arguments. If `dest_name` is not given, the name of the template is
        used. If `dest_dir` is not given, `dest_dir` will be set to the current
        working directory. When `overwrite` is False (default) and the
        destination file already exists, a FileExistsError will be raised.

        Args:
            template_name (str): The name of the template file to be rendered.
            dest_name: The name of the destination file. If not given, the name
                of the template will be used.
            dest_dir: The directory where the destination file will be placed.
                If not given, the current working directory will be used.
            overwrite: If True, the destination file will be overwritten if it
                already exists. If False (default), a FileExistsError will be
                raised if the destination file already exists.
            render_args: Keyword arguments that will be passed to the template
                renderer.

        Raises:
            FileExistsError: If the destination file already exists and
                `overwrite` is False.
        """

        # default name is the same as the template
        if dest_name is None:
            dest_name = pathlib.Path(template_name).name
        dest_name = pathlib.Path(dest_name)

        # default destination directory is the current working directory
        if dest_dir is None:
            dest_dir = pathlib.Path.cwd()
        dest_dir = pathlib.Path(dest_dir)

        dest_path = self.renderer.strip_extension(dest_dir / dest_name)
        if not overwrite and dest_path.exists():
            raise FileExistsError(f"{dest_path} already exists.")

        dest_path.parent.mkdir(parents=True, exist_ok=True)
        template_name = pathlib.Path(template_name).as_posix()
        self.renderer.render(template_name, dest_path,
                             **self._per_template_ctx.get(template_name, {}),
                             **render_args)

    def _validate_destination(self, dest_dir: pathlib.Path):
        """Check if the destination filenames exist.

        Check if the destination filenames exist and raise a FileExistsError if
        any of them already exists.

        Args:
            dest_dir (pathlib.Path): The destination directory.
        """
        for subdir, contents in self.dir_struct.items():
            dest_subdir = dest_dir / subdir
            for file in contents["files"] + contents["templates"]:
                dest_name = self.renderer.strip_extension(dest_subdir / file)
                if dest_name.exists():
                    raise FileExistsError(f"File '{dest_name}' already exists.")


if __name__ == "__main__":
    t = Template("/Users/ssantos/repos/wind-tunnel-demos/lib")
    # pprint.pprint(t.structure)
    kwds = {"flow_velocity": [1, 2, 3], "resolution": 1}
    t.render_dir("treta", overwrite=1, **kwds)

    # t.render_file("templates/constant/triSurface/README.md")
    t.render_file("templates/constant/triSurface/README.md", "tmp.txt")
    t.render_file("templates/constant/triSurface/README.md", "lixo/tmp.txt")
    # t.render_file("templates/constant/triSurface/README.md",
    #               dest_dir="lixo", overwrite=True)
    t.render_file("templates/constant/triSurface/README.md",
                  "subdir",
                  dest_dir="lixo", overwrite=True)
