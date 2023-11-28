"""Resources CLI subcommand."""
import inductiva
from inductiva import _cli


def list_resources(args):
    del args  # unused
    inductiva.resources.machine_groups.list()
    inductiva.resources.storage.get_space_used()


def estimate_machine_cost(args):
    machine_type = args.machine_type
    zone = args.zone
    spot = args.spot

    cost = inductiva.resources.estimate_machine_cost(
        machine_type=machine_type,
        zone=zone,
        spot=spot,
    )
    print(f"Estimated cost of machine: {cost} $/h.")


def register_resources_cli(parser):
    _cli.utils.show_help_msg(parser)
    subparsers = parser.add_subparsers()

    list_subparser = subparsers.add_parser(
        "list", help="List currently active resources")

    cost_subparser = subparsers.add_parser(
        "cost",
        help="Estimate cost of a machine in the cloud",
    )

    cost_subparser.add_argument(
        "machine_type",
        type=str, help="Type of machine to launch")
    cost_subparser.add_argument(
        "-z", "--zone",
        default="europe-west1-b",
        type=str, help="Type of machine to launch")
    cost_subparser.add_argument(
        "--spot",
        default=False,
        type=bool, help="Type of machine to launch")

    # Register function to call when this subcommand is used
    list_subparser.set_defaults(func=list_resources)
    cost_subparser.set_defaults(func=estimate_machine_cost)
