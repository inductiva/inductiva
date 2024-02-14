"""Utils to setup python autocompletion."""
import os
import pathlib
import shutil
from importlib.resources import files

import inductiva
from inductiva import constants


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
    version = inductiva.__version__.replace(".", "-")
    version_dir = constants.LOCAL_LOGGING_DIR / f"v{version}"
    completion_dir = version_dir / "completions"

    if not os.path.exists(version_dir):
        os.makedirs(version_dir)

    if not os.path.exists(completions_dir):
        os.mkdir(completions_dir)

    completion_file = files(
        "inductiva.completions") / f"v{version}" / "_inductiva"

    shutil.copyfile(completion_file, completion_dir "_inductiva")

    lines_to_append = [
        f"fpath=({completion_dir} $fpath)\n",
        "autoload -Uz compinit\n",
        "compinit\n",
    ]

    setup_file_path = version_dir / "inductiva-setup.sh"
    _append_lines_to_file(lines_to_append, setup_file_path)

    lines_to_append = [
        "\n\n# BEGIN: LINES ADDED BY INDUCTIVA THE PACKAGE\n",
        f"source {setup_file_path}\n"
        "# END: LINES ADDED BY THE INDUCTIVA PACKAGE\n\n",
    ]
    zshrc_path = pathlib.Path.home() / ".zshrc"
    _append_lines_to_file(lines_to_append, zshrc_path)


def _append_lines_to_file(lines: List[str], file_path: Union[str,
                                                             pathlib.Path]):
    """
    Append a list of lines to a specified file.

    Parameters:
    - lines (list): List of lines to append to the file.
    - file_path (str or Path): Path to the file.

    Returns:
    - None
    """
    with open(file_path, "a") as f:
        for line in lines:
            f.write(line)
