"""Launch simulation resources via CLI."""

import argparse
from inductiva import resources


def start_machine_group(args):
    """Start a resource."""
    machine_type = args.machine_type
    num_machines = args.num_machines
    data_disk_gb = args.data_disk_gb
    spot = args.spot

    machine = resources.MachineGroup(machine_type=machine_type,
                                     num_machines=num_machines,
                                     data_disk_gb=data_disk_gb,
                                     spot=spot)

    machine.start()
    print(f"{repr(machine)} started.")


def register(parser):
    """Register the launch of resources commands."""

    subparser = parser.add_parser("start",
                                  help="Start computational resources.",
                                  formatter_class=argparse.RawTextHelpFormatter)

    subparser.description = (
        "The `inductiva resources start` command initiates "
        "computational resources.\n"
        "It provides options to specify the machine type, "
        "number of machines,"
        "disk size in GB, and the use of spot instances.\n"
        "This command is suitable"
        "for initiating both individual machines and clusters.\n\n"
        "To terminate resources, use the `inductiva resources"
        " terminate` command.\n")
    subparser.add_argument("machine_type",
                           type=str,
                           help="Machine type to start.")
    subparser.add_argument("-n",
                           "--num_machines",
                           default=1,
                           type=int,
                           help="Number of machines to start.")
    subparser.add_argument("-d",
                           "--data_disk_gb",
                           default=70,
                           type=int,
                           help="Disk size in GB.")
    subparser.add_argument("-s",
                           "--spot",
                           default=True,
                           action="store_true",
                           help="Whether to use spot instances.")

    subparser.set_defaults(func=start_machine_group)
