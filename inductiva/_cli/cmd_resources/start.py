"""Launch simulation resources via CLI."""

from inductiva import resources


def start_machine_group(args):
    """Start a resource."""
    machine_type = args.machine_type
    num_machines = args.num_machines
    disk_size_gb = args.disk_size
    spot = args.spot

    machine = resources.MachineGroup(machine_type=machine_type,
                                     num_machines=num_machines,
                                     disk_size_gb=disk_size_gb,
                                     spot=spot)

    machine.start()
    print(f"{repr(machine)} started.")


def register(parser):
    """Register the launch of resources commands."""

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
    subparser.add_argument("-s",
                           "--spot",
                           default=False,
                           type=bool,
                           help="Whether to use spot instances")

    subparser.set_defaults(func=start_machine_group)
