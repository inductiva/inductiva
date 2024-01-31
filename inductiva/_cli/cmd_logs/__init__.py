"""Register CLI commands for logs."""
import os

from inductiva._cli import loader, utils


def register(root_parser):

    parser = root_parser.add_parser("logs",
                                    help="Stream the logs of a running task.")
    utils.show_help_msg(parser)

    loader.load_commands(parser, os.path.dirname(__file__), package=__name__)
