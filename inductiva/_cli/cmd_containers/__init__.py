"""Register CLI commands for containers."""
import argparse
import os

from .. import loader, utils
from ... import constants


def register(root_parser):
    parser = root_parser.add_parser(
        "containers",
        help="Containers commands",
        formatter_class=argparse.RawTextHelpFormatter)

    parser.description = \
"""
Containers management utilities.

The `inductiva containers` command allows you to manage custom user containers.

It allows users to convert Docker images into Apptainer-compatible .sif files 
and upload them to their Inductiva private storage for use in simulations.
"""
    utils.show_help_msg(parser)

    subparsers = parser.add_subparsers(title="available subcomands")
    loader.load_commands(subparsers,
                         os.path.dirname(__file__),
                         package=__name__,
                         ignores_prefix=constants.LOADER_IGNORE_PREFIX)
