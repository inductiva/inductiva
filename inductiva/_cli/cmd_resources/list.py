"""CLI commands for listing both the active and available machine groups."""

from typing import TextIO
import argparse
import bisect
import sys

from inductiva.resources.machine_types import ProviderType
from inductiva.utils import format_utils
from inductiva import resources, _cli


def pretty_print_machines_info(machines_dict):
    """Format and print the information for given machine dict."""
    emph_formatter = format_utils.get_ansi_formatter()
    header_formatters = [
        lambda x: emph_formatter(x.upper(), format_utils.Emphasis.BOLD)
    ]
    print()
    for family, family_details in machines_dict.items():
        print(f"CPU family: {family}")
        final_table = {"Machine Type": [], "Supported vCPUs": [], "Config": []}
        first_line = True
        for machine_type, details in family_details.items():
            # Used to determine if we write the machine type name
            first_time_type = True

            # Don't want to add an empty line before the first machine type
            if not first_line:
                # Add's an empty line between the machine types
                final_table["Machine Type"].append("")
                final_table["Supported vCPUs"].append("")
                final_table["Config"].append("")
            first_line = False

            for config, vcpus in details.items():
                if first_time_type:
                    final_table["Machine Type"].append(machine_type)
                    # The first time we add the machine type, we don't want to
                    # write its name
                    first_time_type = False
                else:
                    # If we have more than one config for the same machine type
                    # we dont want to repeat the machine type name
                    # (ex c3 standard)
                    final_table["Machine Type"].append("")

                str_vcpus = ", ".join([str(v) for v in vcpus["vcpus"]])
                final_table["Supported vCPUs"].append(str_vcpus)
                final_table["Config"].append(config)

        res_table = format_utils.get_tabular_str(
            final_table,
            header_formatters=header_formatters,
            indentation_level=4)
        print(res_table)

    print("Ex:\tc3d-highcpu-4")
    print("\tc3d-highcpu-8-lssd")


def _replace_empty_lists_with_value(d, value_to_replace):
    """Replace empty lists with a value in a nested dictionary.
    Used to prevent printing of [] in the pretty print.
    """
    for key, value in d.items():
        if isinstance(value, dict):
            _replace_empty_lists_with_value(value,
                                            value_to_replace=value_to_replace)
        elif isinstance(value, list):
            if len(value) == 0:
                d[key] = value_to_replace


def list_machine_types_available(args):
    """List available machine types information per provider and CPU series."""

    provider = args.provider
    machine_family = args.family

    resources_available = resources.machine_types.get_available_machine_types(
        provider, machine_family)

    machines_dict = {}

    for machine in resources_available:
        # Avoid the duplication of adding machine types with spot instances.
        if machine["spot"]:
            continue

        machine_type = machine["machine_type"]
        machine_info = machine_type.split("-")

        family = machine_info[0]
        memory = machine_info[1]
        vcpus = machine_info[2] if len(machine_info) > 2 else None
        config = machine_info[3] if len(machine_info) > 3 else None

        if family not in machines_dict:
            machines_dict[family] = {}
        if memory not in machines_dict[family]:
            machines_dict[family][memory] = {}
        if config not in machines_dict[family][memory]:
            machines_dict[family][memory][config] = {"vcpus": []}

        if vcpus is not None:
            # Sorted insertion of vcpus
            bisect.insort(machines_dict[family][memory][config]["vcpus"],
                          int(vcpus))
    _replace_empty_lists_with_value(machines_dict, "n/a")
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
