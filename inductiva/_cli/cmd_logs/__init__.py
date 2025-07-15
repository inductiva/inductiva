"""Register CLI commands for logs."""
import argparse
import os
import textwrap

from .. import loader, utils
from ... import constants


def register(root_parser):
    parser = root_parser.add_parser(
        "logs",
        help="Stream the logs of a running task.",
        formatter_class=argparse.RawTextHelpFormatter)

    parser.description = textwrap.dedent("""\
        The `inductiva logs` command streams the standard output (STDOUT) of a 
        running task in real-time, useful for monitoring live execution
        progress.
        
        This command will stream the logs being written to STDOUT by a task
        only if the task is still in progress. If the task has finished, or
        transitioned to a status, the logs will not be available for streaming.
        
        Real-time streaming of a running task's standard error (STDERR) is also 
        supported via an argument.
    """)
    utils.show_help_msg(parser)

    loader.load_commands(parser,
                         os.path.dirname(__file__),
                         package=__name__,
                         ignores_prefix=constants.LOADER_IGNORE_PREFIX)
