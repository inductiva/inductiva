"""Get the user's remote storage size in GB via CLI."""

import argparse
from inductiva import storage


def storage_used(unused_args):
    """List the user's remote storage contents."""
    storage.get_space_used()


def register(parser):
    """Register the list user's remote storage command."""

    subparser = parser.add_parser("size",
                                  help="Get the occupied remote storage size.",
                                  formatter_class=argparse.RawTextHelpFormatter)

    subparser.description = (
        "The `inductiva storage size` command calculates "
        "the total size of your data on the platform.\n"
        "It returns the total size in GB of all items in your storage.\n")

    subparser.set_defaults(func=storage_used)
