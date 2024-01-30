"""Mixin for file management.

Provides utility classes and functions for managing files and directories
and rendering template files for generation of simulation inputs.
"""
from typing import Union
import pathlib
import shutil
import glob
import io
import os
import re

from absl import logging
import jinja2

from inductiva.utils import format_utils

TEMPLATE_EXTENSION = ".jinja"
SUFFIX_SEPARATOR = "__"


class FileManager:
    """Class for file management."""
    DEFAULT_ROOT_DIR = "input_dir"
    __root_dir = None  #pylint: disable=invalid-name

    def __check_root_dir(self):
        """Set default root folder if unset."""
        if self.__root_dir is None:
            self.set_root_dir(self.DEFAULT_ROOT_DIR)

    def set_root_dir(self, root_dir=None):
        """Set a root directory for the file manager.
        
        All files managed through the manager will be
        placed inside a newly created folder with the given
        name. If a folder with the same name already exists,
        the newly created one will be appended with the
        "__1" suffix. If folders having the "__N" suffix already
        exist, the given name will be appended with "__{N+1}"
        where N is determined by the directory with the largest N.

        Args:
            root_dir: Path to the root directory.
                If None, an error is raised.
        """

        if root_dir is None:
            raise ValueError("Given root directory cannot be None")

        if os.path.isdir(root_dir):
            if format_utils.getenv_bool(
                    "INDUCTIVA_DISABLE_FILEMANAGER_AUTOSUFFIX", False):
                raise FileExistsError(f"Directory {root_dir} already exists.")
            generated_root_dir = _gen_unique_name(str(root_dir))
            logging.info(
                "Directory %s already exists."
                " Setting root folder to %s.", root_dir, generated_root_dir)
            root_dir = generated_root_dir

        logging.info("Setting root folder to %s.", root_dir)
        os.makedirs(root_dir)
        self.__root_dir = pathlib.Path(root_dir).resolve()

    def get_root_dir(self):
        """Get the active root directory for the file manager."""

        self.__check_root_dir()

        return self.__root_dir

    def add_file(self, source_file, target_file=None, **render_args):
        """Add a file to the root_dir and render in case it is a template file. 

        This method copies the contents of the `source_file` to a `target_file`
        inside the root directory. If the source file is a template file, it is
        rendered using the given `render_args` before copying to the target
        destination.

        If `target_file` is None, the copied file is saved inside the manager's
        root directory with the same name as the source file; if it is a
        template file, the template suffix is removed. 
        A file is considered to be a template file if it ends with the `.jinja`
        suffix. The rendering arguments are ignored if the source file is not a
        template file.

        Args:
            source_file: Path to the source file.
            target_file: Path to the target file.
            render_args: Arguments to render the template file.

        Examples:
        a) Add a file to the root directory "example_a"
            (input_file.txt -> example_a/input_file.txt):
            >>> file_manager = FileManager()
            >>> file_manager.set_root_dir("example_a")
            >>> file_manager.add_file("input_file.txt")
        
        b) Add a file to the root directory "example_b" with a different name
            (input_file.txt -> example_b/renamed_file.txt):
            >>> file_manager = FileManager()
            >>> file_manager.set_root_dir("example_b")
            >>> file_manager.add_file("input_file.txt", "renamed_file.txt")
        
        c) Render a template file to the root directory "example_c" with a
            template argument named "value" (input_file.txt.jinja ->
            example_c/input_file.txt):
            >>> file_manager = FileManager()
            >>> file_manager.set_root_dir("example_c")
            >>> file_manager.add_file("input_file.txt.jinja", value=2)
        """

        self.__check_root_dir()

        if target_file is None:
            new_target_file = pathlib.Path(source_file).name
        else:
            new_target_file = target_file

        new_target_file = os.path.join(self.__root_dir, new_target_file)
        if source_file.endswith(TEMPLATE_EXTENSION):
            if target_file is None:
                new_target_file, _ = os.path.splitext(new_target_file)

            render_file(template_file=source_file,
                        target_file=new_target_file,
                        **render_args)
        else:
            if render_args:
                logging.debug(
                    f"Ignoring render_args since %s"
                    f"isn't a {TEMPLATE_EXTENSION} file", source_file)
            shutil.copy(source_file, new_target_file)

        return new_target_file

    def add_dir(self, source_dir, target_dir=None, **render_args):
        """Add a directory to the root_dir and render all template files.

        This method copies the contents of the `source_dir` to a `target_dir`
        inside the root directory, keeping the structure as the source
        directory. All template files in the source directory are rendered
        using the given `render_args` before copying to the target destination.
        
        If `target_dir` is None, the source directory is copied to the root
        directory. All template files will have the template suffix stripped
        from their name upon copying.
        A file is considered to be a template file if it ends with the `.jinja`
        suffix. The rendering arguments are ignored no template file is
        available in the source directory.

        Args:
            source_dir: Path to the source directory.
            target_dir: Path to the target directory.
            render_args: Arguments to render the template files.

        Examples:
        a) Add a directory to the root directory "example_a"
            (input_dir/ -> example_a/input_dir/):
            >>> file_manager = FileManager()
            >>> file_manager.set_root_dir("example_a")
            >>> file_manager.add_dir("input_dir")

        b) Add directory to the root directory "example_b" with different name
            (input_dir/ -> example_b/renamed_dir/):
            >>> file_manager = FileManager()
            >>> file_manager.set_root_dir("example_b")
            >>> file_manager.add_dir("input_dir", "renamed_dir")
        
        c) Render template directory to the root directory with the template
        arguments named "value" and "name" (input_dir/ -> example_c/input_dir/):
            >>> file_manager = FileManager()
            >>> file_manager.set_root_dir("example_c")
            >>> file_manager.add_dir("input_dir", value=2, name="example")
        """
        self.__check_root_dir()

        if target_dir is None:
            target_dir = "."

        target_dir = os.path.join(self.__root_dir, target_dir)
        shutil.copytree(source_dir,
                        target_dir,
                        dirs_exist_ok=True,
                        symlinks=True)

        if render_args:
            template_files = get_template_files(target_dir)

            for template_file in template_files:
                template_path = os.path.join(target_dir, template_file)
                target_path = template_path.split(TEMPLATE_EXTENSION)[0]

                render_file(template_file=template_path,
                            target_file=target_path,
                            **render_args)
                os.remove(template_path)

        return target_dir


def render_file(template_file, target_file: Union[str, io.StringIO],
                **render_args):
    """Render a file from a template.

    Args:
        template_file: Path to the source file.
        target_file: Path to the target file.
        render_args: Arguments to render the template file.
    """

    source_path = pathlib.Path(template_file)

    source_dir = source_path.parent
    source_file = source_path.name

    environment = jinja2.Environment(loader=jinja2.FileSystemLoader(source_dir))
    template = environment.get_template(source_file)
    stream = template.stream(**render_args)
    stream.dump(target_file)


def get_template_files(src_dir):
    """Get all template files in a directory."""

    template_paths = glob.glob(os.path.join(src_dir, "**",
                                            "*" + TEMPLATE_EXTENSION),
                               recursive=True)

    template_files = [
        os.path.relpath(file_path, start=src_dir)
        for file_path in template_paths
    ]

    return template_files


def _gen_unique_suffix(name, filenames):
    """Generate a suffix for a filename, based on a list of files.
    
    Args:
        name: Name of the file.
        filenames: List of filenames.
    """
    regex = f"^{name}{SUFFIX_SEPARATOR}([0-9]+$)"

    max_suffix = -1
    for filename in filenames:
        if match := re.match(regex, filename):
            max_suffix = max(max_suffix, int(match.group(1)))

    exists_solo = name in filenames
    exists_derived = max_suffix >= 0

    if exists_derived:
        suffix = f"{SUFFIX_SEPARATOR}{max_suffix+1}"
    elif exists_solo:
        suffix = f"{SUFFIX_SEPARATOR}2"
    else:
        suffix = ""
    return suffix


def _gen_unique_name(name):
    """Generates an unique folder name.

    Generates an unique folder name derived from the given name and
    based on the names of the existing folders in the current working directory.
    If no folder exists with the same name as the given name,
    the input name will be returned; otherwise, a suffix is appended
    to the input name to make it unique. 

    See `_gen_unique_suffix` for more details.
    """

    name = name.strip()
    filenames = glob.glob(f"{name}*")
    return name + _gen_unique_suffix(name, filenames)
