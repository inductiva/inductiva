"""Base and specific renderers packages to manage template files."""
from typing import Iterable, Union, IO
from abc import abstractmethod
import pathlib
import os

import jinja2

from inductiva import types, utils


class RendererAdapter:
    """Abstract base class for template renderers."""

    TEMPLATE_EXTENSION = None

    @abstractmethod
    def set_template_dir(self, template_dir: types.PathOrStr) -> None:
        """Set the directory where the templates are located."""
        pass

    @abstractmethod
    def render(self, source_file: Union[types.PathOrStr, IO],
               target_file: types.PathOrStr, **render_args) -> pathlib.Path:
        """Render a file."""
        pass

    @classmethod
    def is_template(cls, file: types.PathOrStr) -> bool:
        """Check if the given file is of template type."""
        if isinstance(file, pathlib.Path):
            return file.suffix == cls.TEMPLATE_EXTENSION
        return file.endswith(cls.TEMPLATE_EXTENSION)

    @classmethod
    def strip_extension(cls, file: types.PathOrStr) -> str:
        """Strip the template extension if the given file is a template."""
        if not cls.is_template(file):
            return file
        if isinstance(file, pathlib.Path):
            return file.with_suffix("")
        return os.path.splitext(file)[0]

    @abstractmethod
    def get_template_files(self) -> Iterable[types.PathOrStr]:
        """Return a list of template files in the given source directory."""
        pass


class JinjaRenderer(RendererAdapter):
    """Jinja2 renderer for .jinja templates."""

    TEMPLATE_EXTENSION = ".jinja"

    def __init__(self, template_dir: types.PathOrStr):
        """Create a renderer adapter for the Jinja2 template engine.
        
        Attributes:
            TEMPLATE_EXTENSION (str): Extension for the template files used by
                Jinja2 templates.
            template_dir (types.PathOrStr): Path to the template directory."""

        self.set_template_dir(template_dir)

    def set_template_dir(self, template_dir: types.PathOrStr) -> None:
        """Set the directory where the templates are located.

        Sets the directory where the templates are located and initializes the
        Jinja2 environment and loader.

        Args:
            template_dir (PathOrStr): The directory where the templates are
                located.
        """
        self.template_dir = pathlib.Path(template_dir)

        loader = jinja2.FileSystemLoader(self.template_dir)
        self.environment = jinja2.Environment(loader=loader,
                                              undefined=jinja2.StrictUndefined)

    def get_template_files(self) -> Iterable[types.PathOrStr]:
        """Return a list of template files in the given source directory."""

        return utils.files.map_dir_files(self.template_dir,
                                         self.TEMPLATE_EXTENSION)

    def render(self, source_file: types.PathOrStr,
               target_file: Union[types.PathOrStr, IO], **render_args):
        """Render a file given a template file and specific render args.
        
        Args:
            target_file (types.PathOrStr): Path to the target file.
            source_file (types.PathOrStr): Path to the source file relative to
                the template_dir.
            render_args: Arguments to render the template file.
        """

        source_path = pathlib.Path(source_file).as_posix()
        template = self.environment.get_template(source_path)
        stream = template.stream(**render_args)

        if isinstance(target_file, pathlib.Path):
            stream.dump(str(target_file))
        else:
            stream.dump(target_file)
