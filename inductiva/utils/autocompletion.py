"""Utils to setup python autocompletion."""
from typing import List, Union, Optional
import os
import pathlib
import shutil

import inductiva
from inductiva import constants

BEGIN_MARKER_LINE = "# INDUCTIVA BEGIN:"
END_MARKER_LINE = "# INDUCTIVA END:"


def setup_zsh_autocompletion():
    """Sets up shell auto-completion for zsh.

    It does this adding to the users .zshrc:

        # BEGIN: LINES ADDED BY INDUCTIVA PACKAGE
        source {setup_file_path}
        # END: LINES ADDED BY THE INDUCTIVA PACKAGE

    That points to a file:
    /.inductiva/v{version}/inductiva-setup.sh. This file than adds the
    autocompletion script to the `fpath` of the user.

    """
    version = inductiva.__version__
    version_dir = constants.LOCAL_LOGGING_DIR / f"v{version}"
    completions_dir = version_dir / "completions"

    os.makedirs(version_dir, exist_ok=True)
    os.makedirs(completions_dir, exist_ok=True)

    completion_file = pathlib.Path(
        inductiva.__path__[0]) / "assets" / "completions" / "zsh" / "_inductiva"

    shutil.copyfile(completion_file, completions_dir / "_inductiva")

    lines_to_append = [
        f"fpath=({completions_dir} $fpath)\n",
        "autoload -Uz compinit\n",
        "compinit\n",
    ]

    setup_file_path = version_dir / "inductiva-setup.sh"
    _append_lines_to_file(lines_to_append, setup_file_path)

    lines_to_append = [
        f"\n\n{BEGIN_MARKER_LINE}\n",
        f"source {setup_file_path}\n"
        f"{END_MARKER_LINE}\n\n",
    ]
    zshrc_path = pathlib.Path.home() / ".zshrc"
    _append_lines_to_file(lines_to_append, zshrc_path, BEGIN_MARKER_LINE,
                          END_MARKER_LINE)


def _append_lines_to_file(lines: List[str],
                          file_path: Union[str, pathlib.Path],
                          remove_between_begin_line: Optional[str] = None,
                          remove_between_end_line: Optional[str] = None):
    """
    Append a list of lines to a specified file.

    Parameters:
    - lines (list): List of lines to append to the file.
    - file_path (str or Path): Path to the file.

    Returns:
    - None
    """
    if os.path.exists(file_path):
        with open(file_path, "r", encoding="utf-8") as f:
            file_content = f.readlines()
    else:
        file_content = []

    if (remove_between_begin_line is not None and
            remove_between_end_line is not None):
        try:
            file_content = _modify_content(file_content,
                                           remove_between_begin_line,
                                           remove_between_end_line)
        except EOFError as e:
            raise EOFError(
                f"Reached end of file while reading {file_path}\n"
                f"Maybe you removed the inductiva marker lines in {file_path}\n"
                f"Please check to see if the lines:\n"
                f"  {BEGIN_MARKER_LINE}\n"
                f"  {END_MARKER_LINE}\n"
                "Exist in your file.") from e

    file_content += lines

    with open(file_path, mode="w", encoding="utf-8") as f:
        for line in file_content:
            f.write(line)


def _modify_content(file_content: List[str], remove_between_begin_line: str,
                    remove_between_end_line: str) -> List[str]:
    """Modify the content of a list of lines by removing lines
    between specified begin and end lines.

    Parameters:
    - file_content (list): List of lines to modify.
    - begin_line (str): The line to start removing content from.
    - end_line (str): The line to end removing content.

    Returns:
    - List of modified lines.

    """
    modified_content = []
    is_inside = False

    for line in file_content:
        if line.startswith(remove_between_begin_line):
            is_inside = True

        if not is_inside:
            modified_content.append(line)

        if line.startswith(remove_between_end_line):
            is_inside = False

    if is_inside:
        raise EOFError

    return modified_content
