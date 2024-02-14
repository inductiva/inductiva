"""Enables autocompletion for the CLI via CLI."""
import argparse

from inductiva.utils.autocompletion import setup_zsh_autocompletion


def enable_auto_complete(args):
    #TODO(augusto): Setup autocompletions for other shells.
    shell = args.shell
    if shell != "zsh":
        raise ValueError("At the moment inductiva supports only autocompletion "
                         "for zsh shells.")
    setup_zsh_autocompletion()
    return 0


def register(parser):
    """Register the list user's tasks command."""

    help_message = "Enables autocomplete for shell commands."
    subparser = parser.add_parser("enable",
                                  help=help_message,
                                  formatter_class=argparse.RawTextHelpFormatter)

    subparser.description = help_message

    group = subparser.add_mutually_exclusive_group()
    group.add_argument("--shell", type=str, help="The shell used. Must be zsh")

    subparser.set_defaults(func=enable_auto_complete)
