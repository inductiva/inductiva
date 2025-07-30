"""Register CLI commands for storage."""
import argparse
import os
import textwrap

from inductiva import constants
from inductiva._cli import loader, utils


def register(root_parser):

    parser = root_parser.add_parser(
        "task-runner",
        help="Task-Runner management utilities.",
        formatter_class=argparse.RawTextHelpFormatter)

    parser.description = textwrap.dedent("""\
        Task-Runner management utilities.
                                         
        A TaskRunner is a local process that manages and executes Inductiva
        tasks. With the Inductiva CLI, you can launch runners locally, so you
        can use your own resources to run your tasks. This allows you to switch
        between local and remote execution, and balance between cost (zero
        compute cost if you are running Tasks on your local machines) and
        performance (access to high-performance machines on the cloud via
        Inductiva). 
        
        The `inductiva task-runner` command allows you to manage your local
        task-runners on the Inductiva platform. It provides utilities including
        launching and terminating task-runners.
    """)

    utils.show_help_msg(parser)

    subparsers = parser.add_subparsers(title="available subcomands")
    loader.load_commands(subparsers,
                         os.path.dirname(__file__),
                         package=__name__,
                         ignores_prefix=constants.LOADER_IGNORE_PREFIX,
                         hides_prefix=constants.LOADER_HIDE_PREFIX)
