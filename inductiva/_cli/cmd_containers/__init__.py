"""Register CLI commands for containers."""
import argparse
import os
import textwrap

from .. import loader, utils
from ... import constants


def register(root_parser):
    parser = root_parser.add_parser(
        "containers",
        help="Containers commands",
        formatter_class=argparse.RawTextHelpFormatter)

    parser.description = textwrap.dedent("""\
        Manage custom simulation containers.

        The `inductiva containers` command provides utilities for managing
        user-defined containers, including converting Docker images to
        Apptainer-compatible  `.sif` files and uploading them to your private
        Inductiva remote storage for use in simulations.
    """)
    utils.show_help_msg(parser)

    subparsers = parser.add_subparsers(title="available subcomands")
    loader.load_commands(subparsers,
                         os.path.dirname(__file__),
                         package=__name__,
                         ignores_prefix=constants.LOADER_IGNORE_PREFIX)
