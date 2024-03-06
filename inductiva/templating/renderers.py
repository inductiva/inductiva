"""Base and specific renderers packages to manage template files."""
from abc import abstractmethod
from typing import Iterable, Union, IO
import pathlib
import shutil
import os

import jinja2

from inductiva import types, utils


class BaseRenderer:
    """Abstract base class for template renderers."""

    TEMPLATE_EXTENSION = None

    @abstractmethod
    def render_dir(self, target_dir: types.Path, source_dir: types.Path,
                   **render_args):
        """Render a template to a file."""
        pass

    @abstractmethod
    def render_file(self, target_file: types.Path,
                    source_file: Union[types.Path, IO], **render_args):
        """Render a file."""
        pass

    @classmethod
    def is_template(cls, file: types.Path) -> bool:
        """Check if the given file is of template type."""
        if isinstance(file, pathlib.Path):
            return file.suffix == cls.TEMPLATE_EXTENSION
        return file.endswith(cls.TEMPLATE_EXTENSION)

    @classmethod
    def strip_extension(cls, file: types.Path) -> str:
        """Strip the template extension if the given file is a template."""
        if not cls.is_template(file):
            return file
        if isinstance(file, pathlib.Path):
            return file.with_suffix("")
        return os.path.splitext(file)[0]

    @classmethod
    def get_template_files(cls, source_dir: types.Path) -> Iterable[types.Path]:
        """Return a list of template files in the given source directory."""

        return utils.files.map_dir_files(source_dir, cls.TEMPLATE_EXTENSION)


class JinjaRenderer(BaseRenderer):
    """Jinja2 renderer for .jinja templates."""

    TEMPLATE_EXTENSION = ".jinja"

    @classmethod
    def render_dir(cls, source_dir: types.Path, target_dir: types.Path,
                   **render_args):
        """Render a directory from a template directory.

        Render the jinja template files within a given source_dir into a 
        target_dir. These files are rendered with the given render_args and
        the rendered files keep the same name without the .jinja extension.

        Note, all non-template files are also copied to the target directory.

        Args:
            source_dir (types.Path): Path to the source directory.
            target_dir (types.Path): Path to the target directory.
            render_args: Arguments to render the template files.
        """
        target_dir = pathlib.Path(target_dir)

        shutil.copytree(source_dir,
                        target_dir,
                        dirs_exist_ok=True,
                        symlinks=True)

        template_dir_struct = cls.get_template_files(target_dir)

        # Render the template files from within the target directory.
        environment = jinja2.Environment(
            loader=jinja2.FileSystemLoader(target_dir),
            undefined=jinja2.StrictUndefined)

        for subdir, contents in template_dir_struct.items():
            for file in contents[cls.TEMPLATE_EXTENSION]:
                target_rel_path = pathlib.Path(subdir) / file
                target_abs_path = target_dir / target_rel_path

                template = environment.get_template(str(target_rel_path))
                stream = template.stream(**render_args)
                stream.dump(str(cls.strip_extension(target_abs_path)))
                target_abs_path.unlink()

    @classmethod
    def render_file(cls, source_file: types.Path,
                    target_file: Union[types.Path, IO], **render_args):
        """Render a file given a template file and specific render args.
        
        Args:
            target_file (types.Path): Path to the target file.
            source_file (types.Path): Path to the source file.
            render_args: Arguments to render the template file.
        """

        source_path = pathlib.Path(source_file)

        environment = jinja2.Environment(loader=jinja2.FileSystemLoader(
            source_path.parent),
                                         undefined=jinja2.StrictUndefined)
        template = environment.get_template(source_path.name)

        stream = template.stream(**render_args)
        if isinstance(target_file, pathlib.Path):
            stream.dump(str(target_file))
        else:
            stream.dump(target_file)
