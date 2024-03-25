"""CLI commands for listing both the active and available machine groups."""

from typing import TextIO
import argparse
import sys

from inductiva import resources, _cli
from inductiva.resources.machine_types import ProviderType


def pretty_print_machines_info(machines_dict):
    """Format and print the information for given machine dict."""

    print()

    for family, memory_types in machines_dict.items():
        print(f"CPU family: {family}")
        for memory, details in memory_types.items():
            vcpus = sorted([int(k) for k in details["vcpus"]], key=int)
            lssd = "(-lssd)" if details["lssd"] else ""
            print(f"   > {family}-{memory}- {(vcpus)} {lssd}")
        print()


def list_machine_types_available(args):
    """List available machine types information per provider and CPU series."""

    provider = args.provider
    machine_family = args.family

    machines_dict = {}

    resources_available = resources.machine_types.get_available_machine_types(
        provider, machine_family)

    for machine in resources_available:
        # Avoid the duplication of adding machine types with spot instances.
        if machine["spot"]:
            continue

        machine_type = machine["machine_type"]
        machine_info = machine_type.split("-")
        machine_family = machine_info[0]
        memory_type = machine_info[1]
        vcpus = machine_info[2] if len(machine_info) > 2 else None
        lssd = True if len(machine_info) > 3 else None

        if machine_family not in machines_dict:
            machines_dict[machine_family] = {}

        if memory_type not in machines_dict[machine_family]:
            machines_dict[machine_family][memory_type] = {
                "vcpus": [],
                "lssd": False
            }

        if vcpus is not None and vcpus not in machines_dict[machine_family][
                memory_type]["vcpus"]:
            machines_dict[machine_family][memory_type]["vcpus"].append(vcpus)

        if lssd is not None:
            machines_dict[machine_family][memory_type]["lssd"] = True

    pretty_print_machines_info(machines_dict)


def list_machine_groups(unused_args, fout: TextIO = sys.stdout):
    """List Resources."""
    resources.machine_groups.list(fout=fout)


def register(parser):
    """Register the list resources commands."""

    subparser = parser.add_parser("available",
                                  help="List available machine types.",
                                  formatter_class=argparse.RawTextHelpFormatter)

    subparser.add_argument("-p",
                           "--provider",
                           type=str,
                           default="gcp",
                           choices=[t.value for t in ProviderType],
                           help="Filter the available types by provider.")
    subparser.add_argument("-f",
                           "--family",
                           default=None,
                           type=str,
                           help="Filter the available types by CPU series.")

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

    _cli.utils.add_watch_argument(subparser)

    subparser.description = (
        "The `inductiva resources list` command provides a snapshot "
        "of your active computational resources.\nFrom machine type "
        "to start time, it gives you a comprehensive overview of "
        "your resources in one place.")
    subparser.set_defaults(func=list_machine_groups)
