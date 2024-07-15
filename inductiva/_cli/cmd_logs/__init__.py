"""Register CLI commands for logs."""
import argparse
import os

from .. import loader, utils
from ... import constants


def register(root_parser):
    parser = root_parser.add_parser(
        "logs",
        help="Stream the logs of a running task.",
        formatter_class=argparse.RawTextHelpFormatter)

    parser.description = (
        "Stream the STDOUT of a running task.\n\n"
        "This command will stream the logs being written to STDOUT by a\n"
        "running task only if the task is still running. If the task has\n"
        "finished, or transitioned to a status other than 'running', the\n"
        "logs will not be available for streaming.\n"
        "To check the status of a task, use:\n"
        "  inductiva tasks list --id <task_id>")

    utils.show_help_msg(parser)

    loader.load_commands(parser,
                         os.path.dirname(__file__),
                         package=__name__,
                         ignores_prefix=constants.LOADER_IGNORE_PREFIX)
