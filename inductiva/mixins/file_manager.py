"""Mixin for file management."""
import glob
import os
import shutil
import pathlib

from absl import logging
import jinja2

from inductiva.utils import files, misc

TEMPLATE_EXTENSION = ".jinja"


class FileManager:
    """Class for file management."""
    DEFAULT_ROOT_DIR = "input_dir"
    __root_dir = None  #pylint: disable=invalid-name

    def __check_root_dir(self):
        if self.__root_dir is None:
            self.set_root_dir(self.DEFAULT_ROOT_DIR)

    def set_root_dir(self, root_dir=None):
        """Set the root directory for the manager."""

        if root_dir is None:
            raise ValueError("Given root directory cannot be None")
        elif os.path.isdir(root_dir):
            generated_root_dir = misc.gen_name(root_dir)
            logging.info(
                "Directory %s already exists."
                " Adding random tag to directory name.", generated_root_dir)
            root_dir = generated_root_dir

        root_dir = files.resolve_path(root_dir)
        os.mkdir(root_dir)
        self.__root_dir = root_dir

    def get_root_dir(self):
        """Get the root directory for the manager."""
        self.__check_root_dir()

        return self.__root_dir

    def add_file(self, source_file, target_file=None, **render_args):
        """Render a file from a template.
        
        If target_file is None, then use the name of the source_file,
        without the template extension.

        Only template files are rendered. Other files are added as is.
        """

        self.__check_root_dir()

        if target_file is None:
            new_target_file = pathlib.Path(source_file).name
        else:
            new_target_file = target_file

        if not source_file.endswith(TEMPLATE_EXTENSION):
            if render_args:
                logging.info(
                    f"Ignoring render_args since %s"
                    f"isn't a {TEMPLATE_EXTENSION} file", source_file)
            shutil.copy(source_file,
                        os.path.join(self.__root_dir, new_target_file))
        else:
            if target_file is None:
                new_target_file = new_target_file.split(TEMPLATE_EXTENSION)[0]

            render_file(source_file=source_file,
                        target_file=os.path.join(self.__root_dir,
                                                 new_target_file),
                        **render_args)

        return new_target_file

    def add_dir(self, source_dir, target_dir=None, **render_args):
        """Render a directory from a template.
        
        Create a new directory inside the root_dir, keeping the same
        structure as the source_dir, and replacing the render_args in
        the template files.

        Only file with the template extension are replaced.
        Other files are added as is.
        """
        self.__check_root_dir()

        if target_dir is None:
            target_dir = "."

        shutil.copytree(source_dir,
                        target_dir,
                        dirs_exist_ok=True,
                        symlinks=True)

        template_files = get_template_files(target_dir)

        for template_file in template_files:
            template_path = os.path.join(target_dir, template_file)
            target_path = template_path.split(TEMPLATE_EXTENSION)[0]

            render_file(source_file=template_path,
                        target_file=target_path,
                        remove_template=True,
                        **render_args)

        return target_dir


def render_file(source_file, target_file, remove_template=False, **render_args):

    source_path = pathlib.Path(source_file)

    source_dir = source_path.parent
    source_file = source_path.name

    environment = jinja2.Environment(loader=jinja2.FileSystemLoader(source_dir))
    template = environment.get_template(source_file)
    stream = template.stream(**render_args)
    stream.dump(target_file)

    if remove_template:
        os.remove(source_path)


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
