import subprocess
import argparse

from inductiva import constants
from inductiva.utils.autocompletion import setup_zsh_autocompletion


def enable_auto_complete(args):
    #TODO(augusto): Setup autocompletions for other shells.
    setup_zsh_autocompletion()
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
