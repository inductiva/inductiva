"""CLI utils."""
import argparse
import pathlib
import shutil
import os
import re

import inductiva
from inductiva import logs

BEGIN_TAG = "# >>> INDUCTIVA BEGIN:"
END_TAG = "# >>> INDUCTIVA END:"


def positive_float(value):
    """Check if the given value is a positive float."""
    try:
        value = float(value)
    except TypeError as _:
        raise argparse.ArgumentTypeError(f"{value} is an invalid type.")
    if value <= 0:
        raise argparse.ArgumentTypeError(f"{value} must be a positive number.")
    return value


def add_watch_argument(subparser):
    """Add watch flag to the subparser that prompts the commands every Nsecs."""
    subparser.set_defaults(watchable=True)
    subparser.add_argument("-w",
                           "--watch",
                           nargs="?",
                           const=2.0,
                           type=positive_float,
                           help="Prompt the command every N seconds.")


def check_running_for_first_time():
    """Checks if it is the first time the cli is being called.

    It does this by checking if there is a folder:

      `$HOME/.inductiva/v{inductiva_version}`

    If there is not returns true and creates that folder.

    """
    version = inductiva.__version__
    dir_name = logs.log.get_logs_file_path().parent / f"v{version}"
    if not os.path.exists(dir_name):
        os.mkdir(dir_name)
        return True
    return False


def user_autocompletion_install_prompt():
    prompt = input("This is the first time running the inductiva cli. "
                   "Would you like to enable autocompletion? (y/[N])?\n"
                   "At the moment only zsh is supported.")
    return prompt.lower() in ["y", "ye", "yes"]


def show_help_msg(parser):
    """Show help message for command if no subcommand is provided."""
    parser.set_defaults(func=lambda _: parser.print_help())


def setup_zsh_autocompletion():
    """Sets up shell auto-completion for zsh.

    It does this by:

    1. Adding the file with completions to:
      `$HOME/.inductiva/v{version}/completions/_inductiva`
    This file already comes with the `inductiva` package.

    2. Adding to the users .zshrc:

        # >>> INDUCTIVA BEGIN:
        # !! Contents within this block are managed by 'inductiva' !!
        # !! Removal might break some functionalities !!
        source {setup_file_path}
        # >>> INDUCTIVA END:

    That points to a file:
    /.inductiva/v{version}/inductiva-setup.sh. This file than adds the
    autocompletion script to the `fpath` of the user.

    """
    version = inductiva.__version__
    version_dir = logs.log.get_logs_file_path().parent / f"v{version}"
    assets_dir = pathlib.Path(inductiva.__path__[0]) / "assets"
    shutil.copytree(assets_dir, version_dir, dirs_exist_ok=True)

    lines_to_append = [
        f"\n\n{BEGIN_TAG}\n",
        "# !! Contents within this block are managed by 'inductiva' !!\n",
        "# !! Removal might break some functionalities !!\n",
        f"source {version_dir / 'inductiva-setup.sh'}\n"
        f"{END_TAG}\n\n",
    ]

    regex = re.compile(f"{BEGIN_TAG}.*{END_TAG}", re.DOTALL)
    zshrc_path = pathlib.Path.home() / ".zshrc"

    with open(zshrc_path, "r", encoding="utf8") as f:
        content = f.read()

    with open(zshrc_path, "w", encoding="utf8") as f:
        content = regex.sub("", content)
        f.seek(0)
        f.truncate()
        lines_as_string = "".join(lines_to_append)
        content += f"\n{lines_as_string}\n"
        f.write(content)
