"""Resources CLI subcommand."""
from inductiva import _cli


def add_start_machines_subparser(parser):
    """Subparser for start command."""

    subparser = parser.add_parser("start", help="Start computational resources")

    subparser.add_argument("machine_type",
                           type=str,
                           help="Machine type to start")
    subparser.add_argument("-n",
                           "--num_machines",
                           default=1,
                           type=int,
                           help="Number of machines to start")
    subparser.add_argument("-d",
                           "--disk_size",
                           default=70,
                           type=int,
                           help="Disk size in GB")
    subparser.add_argument("-z",
                           "--zone",
                           default="europe-west1-b",
                           type=str,
                           help="Zone to start the machines")
    subparser.add_argument("-s",
                           "--spot",
                           default=False,
                           type=bool,
                           help="Whether to use spot instances")

    subparser.set_defaults(func=_cli.machines.start_machine_group)


def add_cost_subparser(parser):
    """Subparser for cost command."""

    subparser = parser.add_parser(
        "cost",
        help="Estimate cost of a machine in the cloud",
    )
    subparser.add_argument("machine_type",
                           type=str,
                           help="Type of machine to launch")
    subparser.add_argument("-z",
                           "--zone",
                           default="europe-west1-b",
                           type=str,
                           help="Type of machine to launch")
    subparser.add_argument("--spot",
                           default=False,
                           type=bool,
                           help="Type of machine to launch")

    subparser.set_defaults(func=_cli.machines.estimate_machine_cost)


def add_terminate_subparser(parser):
    """Subparser for terminate command."""
    subparser = parser.add_parser("terminate", help="Terminate a resource.")
    subparser.add_argument("name",
                           type=str,
                           help="Name of the resource to terminate")
    subparser.set_defaults(func=_cli.machines.terminate_machine_group)


def add_list_subparser(parser):
    """Subparser for list command."""
    subparser = parser.add_parser("list",
                                  help="List currently active resources")
    subparser.set_defaults(func=_cli.machines.list_machine_groups)


def add_available_subparser(parser):
    """Subparser for available command."""
    subparser = parser.add_parser("available",
                                  help="List available machine types")
    subparser.set_defaults(func=_cli.machines.list_machine_types_available)


def register_machines_cli(parser):
    _cli.utils.show_help_msg(parser)
    subparsers = parser.add_subparsers()

    add_list_subparser(subparsers)
    add_available_subparser(subparsers)
    add_start_machines_subparser(subparsers)
    add_cost_subparser(subparsers)
    add_terminate_subparser(subparsers)
