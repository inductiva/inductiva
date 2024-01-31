"""CLI commands for listing both the active and available machine groups."""

from inductiva import resources


def list_machine_types_available(unused_args):
    """List all available machines types."""

    print("Available machine types\n")
    print("machine-type: [cores-available]")
    for (machine_type,
         cores) in resources.machine_types.AVAILABLE_MACHINES.items():
        cores_str = ", ".join(str(core) for core in cores)
        print(f"{machine_type}: [{cores_str}]")

    print("\n E.g. of machine-type: c2-standard-8\n")


def list_machine_groups(unused_args):
    """List Resources."""
    resources.machine_groups.list()


def register(parser):
    """Register the list resources commands."""

    subparser = parser.add_parser("available",
                                  help="List available machine types")
    subparser.set_defaults(func=list_machine_types_available)

    subparser = parser.add_parser("list",
                                  aliases=["ls"],
                                  help="List currently active resources")
    subparser.set_defaults(func=list_machine_groups)
