"""Mixin for file management.

Methods to manage files and directories, within your local computer in
order to prepare them for the simulation process. This mixin can be
used to copy files and directories, and render template files with
the given set of arguments.
"""
import glob
import os
import shutil
import pathlib
import re

from absl import logging
import jinja2

from inductiva.utils import files

TEMPLATE_EXTENSION = ".jinja"


class FileManager:
    """Class for file management."""
    DEFAULT_ROOT_DIR = "input_dir"
    __root_dir = None  #pylint: disable=invalid-name

    def __check_root_dir(self):
        """Check if the root directory is set, and set it otherwise."""
        if self.__root_dir is None:
            self.set_root_dir(self.DEFAULT_ROOT_DIR)

    def set_root_dir(self, root_dir=None):
        """Set a root directory for the file manager.
        
        Args:
            root_dir: Path to the root directory.
                If None, an error is raised.
                If the directory already exists. A new directory
                name is created following a sequence of #1, #2, #3, etc."""

        if root_dir is None:
            raise ValueError("Given root directory cannot be None")
        elif os.path.isdir(root_dir):
            generated_root_dir = gen_name(root_dir)
            logging.info(
                "Directory %s already exists."
                " Adding random tag to directory name.", generated_root_dir)
            root_dir = generated_root_dir

        root_dir = files.resolve_path(root_dir)
        os.mkdir(root_dir)
        self.__root_dir = root_dir

    def get_root_dir(self):
        """Get the root directory for the file manager.
        
        If it hasn't been set yet, a default directory is created."""
        self.__check_root_dir()

        return self.__root_dir

    def add_file(self, source_file, target_file=None, **render_args):
        """Render a file from a template.
        
        If target_file is None, then use the name of the source_file,
        without the template extension. Only template files are rendered. 
        Other files are added as is.

        Args:
            source_file: Path to the source file.
            target_file: Path to the target file.
            render_args: Arguments to render the template file.
        """

        self.__check_root_dir()

        if target_file is None:
            new_target_file = pathlib.Path(source_file).name
        else:
            new_target_file = target_file

        new_target_file = os.path.join(self.__root_dir, new_target_file)
        if not source_file.endswith(TEMPLATE_EXTENSION):
            if render_args:
                logging.info(
                    f"Ignoring render_args since %s"
                    f"isn't a {TEMPLATE_EXTENSION} file", source_file)
            shutil.copy(source_file, new_target_file)
        else:
            if target_file is None:
                new_target_file = new_target_file.split(TEMPLATE_EXTENSION)[0]

            render_file(source_file=source_file,
                        target_file=new_target_file,
                        **render_args)

        return new_target_file

    def add_dir(self, source_dir, target_dir=None, **render_args):
        """Render a directory from a template.
        
        Create a new directory inside the root_dir, keeping the same
        structure as the source_dir, and replacing the render_args in
        the template files. Only file with the template extension are replaced.
        Other files are added as is.

        Args:
            source_dir: Path to the source directory.
            target_dir: Path to the target directory.
            render_args: Arguments to render the template files.
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

                render_file(source_file=template_path,
                            target_file=target_path,
                            **render_args)
                os.remove(template_path)

        return target_dir


def render_file(source_file, target_file, **render_args):
    """Render a file from a template.

    Args:
        source_file: Path to the source file.
        target_file: Path to the target file.
        render_args: Arguments to render the template file.
    """

    source_path = pathlib.Path(source_file)

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


def gen_suffix(name, filenames):
    """Generate a suffix for a filename, based on a list of files.
    
    Args:
        name: Name of the file.
        filenames: List of filenames.
    """
    regex = f"^{name}#([0-9]+$)"

    max_suffix = -1
    for filename in filenames:
        if match := re.match(regex, filename):
            max_suffix = max(max_suffix, int(match.group(1)))

    exists_solo = name in filenames
    exists_derived = max_suffix >= 0

    if exists_derived:
        suffix = f"#{max_suffix+1}"
    elif exists_solo:
        suffix = "#2"
    else:
        suffix = ""
    return suffix


def gen_name(name):
    """Generate a sequential name for a file.
     
    If that name already exists in the current directory, a new one is created
    following a sequential pattern name#1, name#2, etc...
    """

    name = name.strip()
    filenames = glob.glob(f"{name}*")
    return name + gen_suffix(name, filenames)
