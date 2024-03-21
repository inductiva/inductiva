"""Utility class to manage files.

Provides utility classes and functions for managing files and directories
and rendering template files for generation of simulation inputs.
"""
from typing import Iterable, Optional
import logging
import pathlib
import shutil
import glob
import re
import os

from inductiva import types
from inductiva.utils import format_utils

SUFFIX_SEPARATOR = "__"

_logger = logging.getLogger(__name__)


class FileManager:
    """Class for file management."""

    DEFAULT_ROOT_DIR = "input_dir"
    __root_dir = None  # pylint: disable=invalid-name

    def __check_root_dir(self) -> pathlib.Path:
        """Set default root folder if unset."""
        if self.__root_dir is None:
            self.set_root_dir(self.DEFAULT_ROOT_DIR)
        return self.__root_dir

    def set_root_dir(self, root_dir: types.PathOrStr):
        """Set a root directory for the file manager.

        All files managed through the manager will be
        placed inside a newly created folder with the given
        name. If a folder with the same name already exists,
        the newly created one will be appended with the
        "__1" suffix. If folders having the "__N" suffix already
        exist, the given name will be appended with "__{N+1}"
        where N is determined by the directory with the largest N.

        Args:
            root_dir (types.PathOrStr): Path to the root directory.
                If None, an error is raised.
        """

        if root_dir is None:
            raise ValueError("Given root directory cannot be None")
        root_dir = pathlib.Path(root_dir).resolve()

        if root_dir.is_dir():
            if format_utils.getenv_bool(
                    "INDUCTIVA_DISABLE_FILEMANAGER_AUTOSUFFIX", False):
                raise FileExistsError(f"Directory {root_dir} already exists.")
            generated_root_dir = _gen_unique_name(str(root_dir))
            _logger.info(
                "Directory %s already exists."
                " Setting root folder to %s.", root_dir, generated_root_dir)
            root_dir = pathlib.Path(generated_root_dir).resolve()

        _logger.info("Setting root folder to %s.", root_dir)
        root_dir.mkdir(parents=True, exist_ok=True)
        self.__root_dir = root_dir

    def get_root_dir(self):
        """Get the active root directory for the file manager."""

        return self.__check_root_dir()

    def copy_file(self,
                  source_file: types.PathOrStr,
                  target_file: Optional[types.PathOrStr] = None,
                  overwrite: bool = False) -> pathlib.Path:
        """Copy a file to the root_dir.

        This method copies the contents of the `source_file` to a `target_file`
        inside the root directory. If `target_file` is None, the copied file is
        named as the source_file, but not the parent directory structure.

        In case, the file with the same name provided already exists, the
        method by default raises a FileExistsError. To overwrite it anyway, set
        the flag `overwrite` to True.

        Args:
            source_file (types.PathOrStr): Path to the source file.
            target_file (types.PathOrStr - Optional): Path to the target file.
            overwrite (bool): Overwrite the target file if it already exists.
        """

        root_dir = self.__check_root_dir()

        if target_file is None:
            target_file = pathlib.Path(source_file).name
        target_file = root_dir / target_file

        if not overwrite and target_file.exists():
            raise FileExistsError(f"{target_file} already exists.")

        target_file.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy(source_file, target_file)

        return target_file

    def copy_dir(self,
                 source_dir: types.PathOrStr,
                 target_dir: Optional[types.PathOrStr] = None,
                 overwrite: bool = False) -> types.PathOrStr:
        """Copy a directory to the root_dir.

        This method copies the contents of the `source_dir` to a `target_dir`
        inside the root directory, keeping the structure as the source
        directory. If `target_dir` is None, the source directory is copied to
        the root directory.

        If any file in the target directory has the same name as the ones
        about to be copied, the default behavior of this method is to raise a
        FileExistsError. To overwrite the target directory anyway, set the flag
        `overwrite` to True.

        Args:
            source_dir (types.PathOrStr): Path to the source directory.
            target_dir (types.PathOrStr - Optional): Path to the target
                directory.
            overwrite (bool): Overwrite the target directory if  already exists.
        """
        root_dir = self.__check_root_dir()

        if target_dir is None:
            target_dir = "."

        target_dir = root_dir / target_dir
        if not overwrite:
            try:
                self._check_precopy_dir(source_dir, target_dir)
            except FileExistsError as e:
                msg = f"{e}; set `overwrite=True` to overwrite existing files."
                raise FileExistsError(msg) from e

        shutil.copytree(source_dir,
                        target_dir,
                        dirs_exist_ok=True,
                        symlinks=True)

        return target_dir

    def _check_precopy_dir(self, source_dir: pathlib.Path,
                           target_dir: pathlib.Path):
        """Check if any file in source_dir exists in target_dir.

        Check if the target filenames exist and raise a FileExistsError if
        any of them already exists.

        Args:
            source_dir (pathlib.Path): The source directory.
            target_dir (pathlib.Path): The target directory.
        """
        target_dir = pathlib.Path(target_dir)

        for dirpath, _, files in os.walk(source_dir, topdown=True):
            target_subdir = target_dir / os.path.relpath(dirpath, source_dir)
            for file in files:
                if (target_subdir / file).exists():
                    raise FileExistsError(
                        f"File {file} already exists in {target_subdir}")


def _gen_unique_suffix(name: str, filenames: Iterable[str]) -> str:
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


def _gen_unique_name(name: str) -> str:
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
