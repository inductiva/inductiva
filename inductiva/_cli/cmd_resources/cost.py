"""CLI commands to get costs of computational resources."""

from inductiva import resources


def estimate_machine_cost(args):
    """Estimate the cost of a certain machine type."""
    machine_type = args.machine_type
    spot = args.spot
    num_machines = args.num_machines

    cost = resources.estimate_machine_cost(
        machine_type=machine_type,
        spot=spot,
    )

    if num_machines == 1:
        print(f"Estimated cost of machine: {cost} $/h.")
    else:
        print("Estimated total cost (per machine)"
              f": {cost*num_machines} ({cost}) $/h.")


def register(parser):
    """Register the cost estimates commands."""

    subparser = parser.add_parser(
        "cost",
        help="Estimate cost of a machine in the cloud",
    )
    subparser.add_argument("machine_type",
                           type=str,
                           help="Type of machine to launch")
    subparser.add_argument("--spot",
                           default=False,
                           action="store_true",
                           help="Type of machine to launch")

    subparser.add_argument("-n",
                           "--num_machines",
                           default=1,
                           type=int,
                           help="Number of machines to launch")

    subparser.set_defaults(func=estimate_machine_cost)
