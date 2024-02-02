"""Register CLI commands for logs."""
import os

from inductiva._cli import loader, utils
from inductiva import constants


def register(root_parser):

    parser = root_parser.add_parser(
        "logs",
        description="Stream the logs of a running task.",
        help="Stream the logs of a running task.")
    utils.show_help_msg(parser)

    loader.load_commands(parser,
                         os.path.dirname(__file__),
                         package=__name__,
                         ignores_prefix=constants.LOADER_IGNORE_PREFIX)
