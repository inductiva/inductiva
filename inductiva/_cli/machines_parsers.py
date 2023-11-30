"""Resources CLI subcommand."""
from inductiva import _cli


def set_start_machines_subparser(subparsers):
    """Subparser for start command."""

    start_subparser = subparsers.add_parser("start",
                                            help="Start a machine group")

    start_subparser.add_argument("machine_type",
                                 type=str,
                                 help="Machine type to start")
    start_subparser.add_argument("-n",
                                 "--num_machines",
                                 default=1,
                                 type=int,
                                 help="Number of machines to start")
    start_subparser.add_argument("-d",
                                 "--disk_size",
                                 default=70,
                                 type=int,
                                 help="Disk size in GB")
    start_subparser.add_argument("-z",
                                 "--zone",
                                 default="europe-west1-b",
                                 type=str,
                                 help="Zone to start the machines")
    start_subparser.add_argument("-s",
                                 "--spot",
                                 default=False,
                                 type=bool,
                                 help="Whether to use spot instances")

    start_subparser.set_defaults(func=_cli.machines.start_machine_group)


def set_cost_subparser(subparsers):
    """Subparser for cost command."""

    cost_subparser = subparsers.add_parser(
        "cost",
        help="Estimate cost of a machine in the cloud",
    )
    cost_subparser.add_argument("machine_type",
                                type=str,
                                help="Type of machine to launch")
    cost_subparser.add_argument("-z",
                                "--zone",
                                default="europe-west1-b",
                                type=str,
                                help="Type of machine to launch")
    cost_subparser.add_argument("--spot",
                                default=False,
                                type=bool,
                                help="Type of machine to launch")

    cost_subparser.set_defaults(func=_cli.machines.estimate_machine_cost)


def set_terminate_subparser(subparsers):
    """Subparser for terminate command."""
    terminate_subparser = subparsers.add_parser(
        "terminate", help="Terminate a machine-group")
    terminate_subparser.add_argument(
        "name", type=str, help="Name of the machine group to terminate")
    terminate_subparser.set_defaults(func=_cli.machines.terminate_machine_group)


def set_list_subparser(subparsers):
    """Subparser for list command."""
    list_subparser = subparsers.add_parser(
        "list", help="List currently active resources")
    list_subparser.set_defaults(func=_cli.machines.list_machine_groups)


def set_available_subparser(subparsers):
    """Subparser for available command."""
    available_subparser = subparsers.add_parser(
        "available", help="List available machine types")
    available_subparser.set_defaults(
        func=_cli.machines.list_machine_types_available)


def register_machines_cli(parser):
    _cli.utils.show_help_msg(parser)
    subparsers = parser.add_subparsers()

    set_list_subparser(subparsers)
    set_available_subparser(subparsers)
    set_start_machines_subparser(subparsers)
    set_cost_subparser(subparsers)
    set_terminate_subparser(subparsers)
