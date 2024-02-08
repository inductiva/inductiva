"""CLI commands for listing both the active and available machine groups."""

import argparse

from inductiva import resources


def list_machine_types_available(unused_args):
    """List all available machines types."""

    resources_available = resources.machine_types.get_available()

    # Fetch resources for each provider
    print()
    for _, provider_resources in resources_available.items():
        description = provider_resources["description"]
        print(f"{description}\n")
        # Fetch the available machine CPU series
        for cpu_series, series_info in provider_resources["cpu-series"].items():
            description = series_info["description"]
            print(f"{cpu_series}: {description}")
            # Fetch the available RAM types and vCPUs info
            for ram_type, type_info in series_info["types"].items():
                _ = type_info["description"]
                vcpus = type_info["vcpus"]
                print(f"  > {cpu_series}-{ram_type}: {vcpus}")
            print()


def list_machine_groups(unused_args):
    """List Resources."""
    resources.machine_groups.list()


def register(parser):
    """Register the list resources commands."""

    subparser = parser.add_parser("available",
                                  help="List available machine types.",
                                  formatter_class=argparse.RawTextHelpFormatter)

    subparser.description = (
        "The `inductiva available` command provides a utility for listing all\n"
        "available machine types along with the number of cores they have.\n\n")

    subparser.set_defaults(func=list_machine_types_available)

    subparser = parser.add_parser("list",
                                  aliases=["ls"],
                                  help="List currently active resources.",
                                  formatter_class=argparse.RawTextHelpFormatter)

    subparser.description = (
        "The `inductiva resources list` command provides a snapshot "
        "of your active computational resources.\nFrom machine type "
        "to start time, it gives you a comprehensive overview of "
        "your resources in one place.")
    subparser.set_defaults(func=list_machine_groups)
