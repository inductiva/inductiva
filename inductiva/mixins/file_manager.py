"""Mixin for file management."""
import os
import shutil
import pathlib
from absl import logging

from inductiva.utils import files, misc, templates


class FileManager:
    """Class for file management."""
    __root_dir = None  #pylint: disable=invalid-name

    def __check_root_dir(self):
        if self.__root_dir is None:
            raise ValueError("Root directory not set.")

    def set_root_dir(self, root_dir=None):
        """Set the root directory for the manager."""

        if root_dir is None:
            logging.info("Root directory not set. "
                         "Setting default root_dir to be `input_dir`.")
            root_dir = files.resolve_path(
                f"input_dir-{misc.create_random_tag()}")
        elif os.path.isdir(root_dir):
            logging.info(
                "Directory %s already exists."
                " Adding random tag to directory name.", root_dir)
            root_dir = files.resolve_path(
                f"{root_dir}-{misc.create_random_tag(size=5)}")

        os.mkdir(root_dir)
        self.__root_dir = root_dir

    def get_root_dir(self):
        """Get the root directory for the manager."""
        self.__check_root_dir()

        return self.__root_dir

    def add_file(self, source_file, target_file=None, **render_args):
        """Render a file from a template.
        
        If target_file is None, then use the name of the source_file,
        without the .jinja extension.

        Only file.jinja are rendered. Other files are added as is.
        """

        self.__check_root_dir()

        if target_file is None:
            new_target_file = pathlib.Path(source_file).name
        else:
            new_target_file = target_file

        if not source_file.endswith(".jinja"):
            if render_args:
                logging.info(
                    "Ignoring render_args since %s"
                    "isn't a jinja file", source_file)
            shutil.copy(source_file,
                        os.path.join(self.__root_dir, new_target_file))
        else:
            if target_file is None:
                new_target_file = new_target_file.split(".jinja")[0]

            templates.render_file(source_file=source_file,
                                  target_file=os.path.join(
                                      self.__root_dir, new_target_file),
                                  **render_args)

        return new_target_file

    def add_dir(self, source_dir, target_dir=None, **render_args):
        """Render a directory from a template.
        
        Create a new directory inside the root_dir, keeping the same
        structure as the source_dir, and replacing the render_args in
        the jinja files.

        Only file.jinja are replaced. Other files are added as is.
        """
        self.__check_root_dir()

        if target_dir is None:
            target_dir = "."

        templates.render_directory(source_dir=source_dir,
                                   target_dir=os.path.join(
                                       self.__root_dir, target_dir),
                                   **render_args)
