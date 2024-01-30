"""CLI commands to terminate computational resources."""

from inductiva import resources


def terminate_machine_group(args):
    """Terminate a machine group from a given name."""
    machine_name = args.name

    print("Terminating MachineGroup... If exists.")
    machines_list = resources.machine_groups.get()

    for machine in machines_list:
        if machine.name == machine_name:
            machine.terminate()
            print(f"Terminated MachineGroup: {machine_name}")
            return

    print(f"MachineGroup {machine_name} not found.")


def register(parser):
    """Register the terminate command for the resources."""

    subparser = parser.add_parser("terminate", help="Terminate a resource.")
    subparser.add_argument("name",
                           type=str,
                           help="Name of the resource to terminate")
    subparser.set_defaults(func=terminate_machine_group)
