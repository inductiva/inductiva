"""CLI commands for listing both the active and available machine groups."""

from typing import TextIO
import argparse
import sys

from inductiva.resources.utils import ProviderType
from inductiva.resources import machine_groups
from inductiva.utils import format_utils
from inductiva import resources, _cli


def pretty_print_machines_info(machines_dict):
    """Format and print the information for given machine dict."""
    emph_formatter = format_utils.get_ansi_formatter()
    header_formatters = [
        lambda x: emph_formatter(x.upper(), format_utils.Emphasis.BOLD)
    ]

    final_table = {
        "Machine Type": [],
        "vCPUS": [],
        "GPUS": [],
        "Memory (GB)": [],
        "Price/Hour (USD)": [],
        "Zone": []
    }
    first_line = True
    for machine_type, zones in machines_dict.items():
        if not first_line:
            # Add's an empty line between the machine types
            final_table["Machine Type"].append("")
            final_table["vCPUS"].append("")
            final_table["GPUS"].append("")
            final_table["Memory (GB)"].append("")
            final_table["Price/Hour (USD)"].append("")
            final_table["Zone"].append("")
        first_line = False

        for zone, details in zones.items():
            for i, vcpu in enumerate(details["vcpus"]):
                final_table["Machine Type"].append(machine_type)
                final_table["vCPUS"].append(vcpu)
                final_table["GPUS"].append(
                    f"{details['gpus'][i]} x {details['gpu_name']}"
                    if details["gpus"] else "n/a")
                final_table["Memory (GB)"].append(details["memory"])
                final_table["Price/Hour (USD)"].append(details["price"])
                final_table["Zone"].append(zone)

    res_table = format_utils.get_tabular_str(
        final_table, header_formatters=header_formatters, indentation_level=4)
    print(res_table)


def list_machine_types_available(args):
    """List available machine types information per provider and CPU series."""

    provider = args.provider
    machine_family = args.family

    resources_available = resources.get_available_machine_types(
        provider, machine_family)
    resources_available.sort(
        key=lambda x: (x.machine_type.split("-")[:-1], x.num_vcpus))

    machines_dict = {}

    for machine in resources_available:
        # Avoid the duplication of adding machine types with spot instances.
        if machine["spot"]:
            continue

        machine_type = machine["machine_type"]

        if machine_type not in machines_dict:
            machines_dict[machine_type] = {}

        memory = machine.ram_gb
        price = machine.price
        vcpus = machine.num_vcpus
        gpus = machine.num_gpus if machine.num_gpus else None
        gpu_name = machine.gpu_name if machine.num_gpus else None
        zone = machine.zone

        machines_dict[machine_type][zone] = {
            "vcpus": [],
            "gpus": [],
            "memory": memory,
            "price": price,
            "gpu_name": gpu_name
        }

        machines_dict[machine_type][zone]["vcpus"].append(int(vcpus))
        if gpus is not None:
            machines_dict[machine_type][zone]["gpus"].append(int(gpus))
    pretty_print_machines_info(machines_dict)


def _machine_group_list_to_str(machine_group_list) -> str:
    """Returns a string representation of a list of machine groups."""
    columns = [
        "Name", "Machine Type", "Elastic", "Type", "# machines",
        "Data Size in GB", "Spot", "Created at (UTC)", "Idle Time",
        "Max Cost ($/hour)"
    ]
    rows = []

    for machine_group in machine_group_list:
        is_elastic = False
        resource_type = machine_groups.ResourceType.STANDARD.value
        spot = machine_group.spot if hasattr(machine_group, "spot") else False

        if isinstance(machine_group, resources.ElasticMachineGroup):
            is_elastic = True
        else:
            if isinstance(machine_group, resources.MPICluster):
                resource_type = machine_groups.ResourceType.MPI.value
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
        "Created at (UTC)": [format_utils.datetime_formatter],
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

    machine_group_list = resources.get()
    if len(machine_group_list) != 0:
        print("Active Resources:", file=fout)
        print(_machine_group_list_to_str(machine_group_list), file=fout, end="")
    else:
        print("No active computational resources found.", file=fout)


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
