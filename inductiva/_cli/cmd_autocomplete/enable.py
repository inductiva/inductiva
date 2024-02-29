"""Enables autocompletion for the CLI via CLI."""
import argparse

from inductiva import _cli


def enable_auto_complete(args):
    #TODO(augusto): Setup autocompletions for other shells.
    shell = args.shell
    if shell != "zsh":
        raise ValueError("At the moment inductiva supports only autocompletion "
                         "for zsh shells.")
    _cli.utils.setup_zsh_autocompletion()
    return 0


def register(parser):
    """Enable autocompletion for the cli."""

    help_message = "Enables autocomplete for shell commands."
    subparser = parser.add_parser("enable",
                                  help=help_message,
                                  formatter_class=argparse.RawTextHelpFormatter)

    subparser.description = help_message

    group = subparser.add_mutually_exclusive_group()
    group.add_argument("--shell",
                       type=str,
                       choices=["zsh"],
                       help="The shell used. Must be zsh")

    subparser.set_defaults(func=enable_auto_complete)
