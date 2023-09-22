"""Resources CLI subcommand."""
import inductiva
from inductiva import _cli


def list_resources(args):
    del args  # unused
    inductiva.resources.machine_groups.list()
    inductiva.resources.storage.get_space_used()


def register_resources_cli(parser):
    _cli.utils.show_help_msg(parser)
    subparsers = parser.add_subparsers()

    list_subparser = subparsers.add_parser(
        "list", help="List currently active resources")

    # Register function to call when this subcommand is used
    list_subparser.set_defaults(func=list_resources)
