"""Utils to setup python autocompletion."""
import os
import pathlib
import subprocess

import inductiva
from inductiva import constants


def setup_zsh_autocompletion():
    version_dir = constants.LOCAL_LOGGING_DIR / f"v{inductiva.__version__}"
    completion_file = version_dir / "completions" / "_inductiva"
    os.makedirs(completion_file.parent, exist_ok=True)
    subprocess.run(f"inductiva --print-completion zsh | tee {completion_file}",
                   shell=True,
                   stdout=subprocess.DEVNULL,
                   stderr=subprocess.DEVNULL)

    lines_to_append = [
        f"fpath=({completion_file.parent} $fpath)\n",
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
