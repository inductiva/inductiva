"""CLI commands to get costs of computational resources."""

import argparse

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

    total_cost = cost * num_machines

    if num_machines == 1:
        print(f"Estimated cost of machine: {cost:.3f} $/h.")
    else:
        print("Estimated total cost (per machine)"
              f": {total_cost:.3f} ({cost:.3f}) $/h.")


def register(parser):
    """Register the cost estimates commands."""

    subparser = parser.add_parser(
        "cost",
        help="Estimate cost of a machine in the cloud.",
        formatter_class=argparse.RawTextHelpFormatter)

    subparser.description = (
        "The `inductiva cost` command provides a utility"
        " for estimating the cost\n"
        "of a machine in the cloud. It allows you to specify"
        " the type of machine\n"
        "and the number of machines, and it can calculate the"
        " cost for spot instances.\n\n")
    subparser.add_argument("machine_type",
                           type=str,
                           help="Type of machine to launch.")
    subparser.add_argument("--spot",
                           default=False,
                           action="store_true",
                           help="Type of machine to launch.")

    subparser.add_argument("-n",
                           "--num_machines",
                           default=1,
                           type=int,
                           help="Number of machines to launch.")

    subparser.set_defaults(func=estimate_machine_cost)
