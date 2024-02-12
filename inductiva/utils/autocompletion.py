"""Utils to setup python autocompletion."""
import os
import pathlib
import shutil
from importlib.resources import files

import inductiva
from inductiva import constants


def setup_zsh_autocompletion():
    version = inductiva.__version__.replace(".", "-")
    version_dir = constants.LOCAL_LOGGING_DIR / f"v{version}"
    completion_dir = version_dir / "completions"
    os.makedirs(completion_dir.parent, exist_ok=True)

    dir_with_completions = files("inductiva.completions") / f"v{version}"
    if not os.path.exists(version_dir / "completions"):
        os.mkdir(version_dir / "completions")
    shutil.copyfile(dir_with_completions / "_inductiva",
                    version_dir / "completions" / "_inductiva")

    lines_to_append = [
        f"fpath=({completion_dir} $fpath)\n",
        "autoload -Uz compinit\n",
        "compinit\n",
    ]

    setup_file_path = version_dir / "inductiva-setup.sh"
    with open(setup_file_path, "a") as f:
        for line in lines_to_append:
            f.write(line)

    lines_to_append = [
        "\n\n# BEGIN: LINES ADDED BY INDUCTIVA PACKAGE\n",
        f"source {setup_file_path}\n"
        "# END: LINES ADDED BY THE INDUCTIVA PACKAGE\n\n",
    ]
    zshrc_path = pathlib.Path.home() / ".zshrc"
    with open(zshrc_path, "a") as f:
        for line in lines_to_append:
            f.write(line)
