"""Get the user's remote storage size in GB via CLI."""

from inductiva import storage

def storage_used(args):
    """List the user's remote storage contents."""
    del args
    storage.get_space_used()

def register(parser):
    """Register the list user's remote storage command."""

    subparser = parser.add_parser("size",
        help="Get the occupied remote storage size.")

    subparser.set_defaults(func=storage_used)
