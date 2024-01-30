"""Register CLI commands for storage."""
import os

from inductiva._cli import loader, utils


def register(root_parser):

    parser = root_parser.add_parser("tasks",
                help="Manage your tasks.")
    utils.show_help_msg(parser)

    subparsers = parser.add_subparsers()
    loader.load_commands(subparsers,
                         os.path.dirname(__file__),
                         package=__name__)
