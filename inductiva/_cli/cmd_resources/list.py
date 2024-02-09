"""CLI commands for listing both the active and available machine groups."""

import argparse

from inductiva import resources


def print_cpu_series_info(cpu_series: str,
                          cpu_series_info: dict,
                          verbose: bool = False):
    """Format and print the information for given CPU series."""

    cpu_series_description = cpu_series_info["description"]
    print(f"{cpu_series}: {cpu_series_description}")
    for ram_type, type_info in cpu_series_info["types"].items():
        description = ("-> " + type_info["description"]) if verbose else ""
        vcpus = str(type_info["vcpus"])
        cpu_type = f"{cpu_series}-{ram_type}-"
        print(f"  > {cpu_type:13} {vcpus:45} {description}")
    print()


def list_machine_types_available(args):
    """List available machine types information per provider and CPU series."""

    provider = args.provider
    description = args.description
    cpu_series = args.series

    resources_available = resources.machine_types.get_available()
    provider_resources = resources_available[provider]
    print(provider_resources["description"] + "\n")

    if cpu_series is not None:
        cpu_series_info = provider_resources["cpu-series"][cpu_series]
        print_cpu_series_info(cpu_series, cpu_series_info, description)
    else:
        # Print all available machine CPU series information
        for cpu_series, series_info in provider_resources["cpu-series"].items():
            print_cpu_series_info(cpu_series, series_info, description)


def list_machine_groups(unused_args):
    """List Resources."""
    resources.machine_groups.list()


def register(parser):
    """Register the list resources commands."""

    subparser = parser.add_parser("available",
                                  help="List available machine types.",
                                  formatter_class=argparse.RawTextHelpFormatter)

    subparser.add_argument("-p",
                           "--provider",
                           type=str,
                           default="gcp",
                           choices=["gcp"],
                           help="Filter the available types by provider.")
    subparser.add_argument("-s",
                           "--series",
                           default=None,
                           type=str,
                           help="Filter the available types by CPU series.")
    subparser.add_argument("-d",
                           "--description",
                           action="store_true",
                           help="Show a description of the machine types.")

    subparser.description = (
        "The `inductiva available` command provides a utility for listing the\n"
        "available machine types by provider (default: gcp) and CPU series.\n"
        "The list includes a description of the memory types and vCPUs "
        "available.\n\n")

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
