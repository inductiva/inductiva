"""Resources CLI subcommand."""
import inductiva


def list_resources(args):
    del args  # unused
    inductiva.resources.machine_groups.list()
    inductiva.resources.storage.get_space_used()


def register_resources_cli(parser):
    subparsers = parser.add_subparsers()
    list_subparser = subparsers.add_parser(
        "list", help="List currently active resources.")

    # Register function to call when this subcommand is used
    list_subparser.set_defaults(func=list_resources)
