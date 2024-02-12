import subprocess
import argparse

from inductiva import constants
from inductiva.utils.autocompletion import setup_zsh_autocompletion


def enable_auto_complete(args):
    print(f"The shell called is {args.shell}")
    shell = args.shell
    if shell not in ["zsh", "bash"]:
        raise ValueError("The provided shell must be either zsh or bash.")

    if shell == "zsh":
        setup_zsh_autocompletion()
        # completion_dir = constants.LOCAL_LOGGING_DIR / f"v{inductiva.__version__}" / "completions" / "_inductiva"
        # os.makedirs(completion_dir.parent, exist_ok=True)
        # subprocess.run(
        #     f"inductiva --print-completion zsh | tee {completion_file}",
        #     shell=True)

        # lines_to_append = [
        #     f"export FPATH=$FPATH:$HOME/{completion_file.parent}\n",
        #     "autoload -Uz compinit", "compinit"
        # ]
        # zshrc_path = Path.home() / ".zshrc"
        # with open(zshrc_path, "a") as file:
        #     for line in lines_to_append:
        #         file.write(line)

    return 0


def register(parser):
    """Register the list user's tasks command."""

    subparser = parser.add_parser("enable",
                                  help="Enables autocomplete",
                                  formatter_class=argparse.RawTextHelpFormatter)

    subparser.description = ("Enables sutocomplete for shell commands.")

    group = subparser.add_mutually_exclusive_group()
    group.add_argument("--shell",
                       type=str,
                       help="The shell used, eitheir bash or zsh")

    subparser.set_defaults(func=enable_auto_complete)
