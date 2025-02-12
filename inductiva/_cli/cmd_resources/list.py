"""CLI commands for listing both the active and available machine groups."""

from typing import TextIO
import argparse
import bisect
import sys

from inductiva.resources.machine_types import ProviderType
from inductiva.resources import machines_base
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
        if family not in ["a3", "g2"]:
            print(f"CPU family: {family}")
        else:
            print(f"GPU family: {family}")
        final_table = {
            "Machine Type": [],
            "Suffix": [],
            "Supported vCPUs": [],
            "Supported GPUs": [],
            "Config": []
        }
        first_line = True
        for machine_type, details in family_details.items():
            # Don't want to add an empty line before the first machine type
            if not first_line:
                # Add's an empty line between the machine types
                final_table["Machine Type"].append("")
                final_table["Suffix"].append("")
                final_table["Supported vCPUs"].append("")
                final_table["Supported GPUs"].append("")
                final_table["Config"].append("")
            first_line = False

            first_entry = True
            for config, info in details.items():
                if family == "a3":
                    for i, gpu in enumerate(info["gpus"]):
                        final_table["Machine Type"].append(
                            machine_type if first_entry else "")
                        final_table["Suffix"].append(f"{gpu}g")
                        final_table["Supported vCPUs"].append(info["vcpus"][i])
                        final_table["Supported GPUs"].append(gpu)
                        final_table["Config"].append(config)
                        first_entry = False
                else:
                    for i, vcpu in enumerate(info["vcpus"]):
                        final_table["Machine Type"].append(
                            machine_type if first_entry else "")
                        final_table["Suffix"].append(vcpu)
                        final_table["Supported vCPUs"].append(vcpu)
                        if info["gpus"]:
                            final_table["Supported GPUs"].append(
                                info["gpus"][i])
                        else:
                            final_table["Supported GPUs"].append("n/a")
                        final_table["Config"].append(config)
                        first_entry = False

        res_table = format_utils.get_tabular_str(
            final_table,
            header_formatters=header_formatters,
            indentation_level=4)
        print(res_table)

    print("Pattern: family-machine type-suffix-config")
    print("Ex:\t c3d-highcpu-4")
    print("\t c3d-highcpu-8-lssd")
    print("\t a3-highgpu-1g")


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
        vcpus = machine.num_cpus
        gpus = machine.num_gpus if machine.num_gpus else None
        config = machine_info[3] if len(machine_info) > 3 else None

        if family not in machines_dict:
            machines_dict[family] = {}
        if memory not in machines_dict[family]:
            machines_dict[family][memory] = {}
        if config not in machines_dict[family][memory]:
            machines_dict[family][memory][config] = {
                "vcpus": [],
                "gpus": [],
            }

        if vcpus is not None:
            # Sorted insertion of vcpus
            bisect.insort(machines_dict[family][memory][config]["vcpus"],
                          int(vcpus))
        if gpus is not None:
            # Sorted insertion of vcpus
            bisect.insort(machines_dict[family][memory][config]["gpus"],
                          int(gpus))
    pretty_print_machines_info(machines_dict)


def _machine_group_list_to_str(machine_group_list) -> str:
    """Returns a string representation of a list of machine groups."""
    columns = [
        "Name", "Machine Type", "Elastic", "Type", "# machines",
        "Data Size in GB", "Spot", "Started at (UTC)", "Idle Time",
        "Max Cost ($/hour)"
    ]
    rows = []

    for machine_group in machine_group_list:
        is_elastic = False
        resource_type = machines_base.ResourceType.STANDARD.value
        spot = machine_group.spot if hasattr(machine_group, "spot") else False

        if isinstance(machine_group, resources.ElasticMachineGroup):
            is_elastic = True
        else:
            if isinstance(machine_group, resources.MPICluster):
                resource_type = machines_base.ResourceType.MPI.value
        num_active_machines = machine_group.active_machines_to_str()

        idle_time = f"{machine_group.idle_time}"
        if machine_group.max_idle_time:
            idle_time += f"/{machine_group.max_idle_time}"

        rows.append([
            machine_group.name, machine_group.machine_type, is_elastic,
            resource_type, num_active_machines, machine_group.data_disk_gb,
            spot, machine_group.create_time, idle_time,
            machine_group.quota_usage.get("max_price_hour")
        ])

    formatters = {
        "Started at (UTC)": [format_utils.datetime_formatter],
    }

    emph_formatter = format_utils.get_ansi_formatter()
    header_formatters = [
        lambda x: emph_formatter(x.upper(), format_utils.Emphasis.BOLD)
    ]

    return format_utils.get_tabular_str(rows,
                                        columns,
                                        formatters=formatters,
                                        header_formatters=header_formatters)


def list_machine_groups(_, fout: TextIO = sys.stdout):
    # pylint: disable=line-too-long
    """Lists all active resources info.

    This method queries Google Cloud for active resources
    that belong to the user. It outputs all the information relative to
    each resource as folllows:

    Active Resources:
                                        Name         VM Type    Elastic   # machines    Disk Size in GB       Spot   Started at (UTC)
    api-6359c03d-c4f9-479f-8b11-ba1f8f55a58c   e2-standard-4      False            3                 40      False   10 Oct, 13:40:50
    api-db2046cf-a6fc-4124-926c-1a24329da5ea   e2-standard-4       True          2/4                 40      False   10 Oct, 12:43:03
    """
    # pylint: enable=line-too-long

    machine_group_list = resources.machine_groups.get()
    if len(machine_group_list) != 0:
        print("Active Resources:", file=fout)
        print(_machine_group_list_to_str(machine_group_list), file=fout, end="")
    else:
        print("No active computational resources found.", file=fout, end="")


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
                           nargs="+",
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
