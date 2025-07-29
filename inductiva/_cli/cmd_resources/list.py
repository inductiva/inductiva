"""CLI commands for listing both the active and available machine groups."""

import textwrap
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

    subparser.description = textwrap.dedent("""\
        The `inductiva resources available` command provides a utility for
        listing the available machine types by provider (default: `GCP`) and
        CPU series. The list includes a description of the memory types and
        vCPUs available.
    """)

    subparser.epilog = textwrap.dedent("""\
        examples:
            # List all available machine types
            $ inductiva resources available

            MACHINE TYPE            VCPUS     GPUS                     MEMORY (GB)     PRICE/HOUR (USD)     ZONE
            a2-highgpu-1g           12        1 x NVIDIA A100 (40gb)   85.0            2.353126             us-central1-a
            a2-highgpu-1g           12        1 x NVIDIA A100 (40gb)   85.0            3.747972             europe-west4-b
            
            a2-highgpu-2g           24        2 x NVIDIA A100 (40gb)   170.0           4.706252             us-central1-a
            a2-highgpu-2g           24        2 x NVIDIA A100 (40gb)   170.0           7.495944             europe-west4-b
            
            a2-highgpu-4g           48        4 x NVIDIA A100 (40gb)   340.0           9.412504             us-central1-a
            a2-highgpu-4g           48        4 x NVIDIA A100 (40gb)   340.0           14.991888            europe-west4-b
            
            a2-highgpu-8g           96        8 x NVIDIA A100 (40gb)   680.0           18.825008            us-central1-a
            a2-highgpu-8g           96        8 x NVIDIA A100 (40gb)   680.0           29.983776            europe-west4-b
            ...
            
            # List the available machine types of the c3d family
            $ inductiva resources available -f c3d

            MACHINE TYPE            VCPUS     GPUS     MEMORY (GB)     PRICE/HOUR (USD)     ZONE
            c3d-highcpu-4           4         n/a      8.0             0.16508132           europe-west1-b
            
            c3d-highcpu-8           8         n/a      16.0            0.33016264           europe-west1-b
            
            c3d-highcpu-16          16        n/a      32.0            0.66032528           europe-west1-b
            
            c3d-highcpu-30          30        n/a      59.0            1.233750645          europe-west1-b
            
            c3d-highcpu-60          60        n/a      118.0           2.46750129           europe-west1-b
            
            c3d-highcpu-90          90        n/a      177.0           3.701251935          europe-west1-b
            
            c3d-highcpu-180         180       n/a      354.0           7.40250387           europe-west1-b
            
            c3d-highcpu-360         360       n/a      708.0           14.80500774          europe-west1-b
            ...
    """)

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

    subparser.epilog = textwrap.dedent("""\
        examples:
            $ inductiva resources list
            Active Resources:

            NAME                            MACHINE TYPE     ELASTIC     TYPE       # MACHINES     DATA SIZE IN GB     SPOT     CREATED AT (UTC)     IDLE TIME      MAX COST ($/HOUR)
            api-tgowxa5pdqxoz3kqtzdyuxay2   c2-standard-4    False       standard   0/1            10                  True     24/07, 14:37:00      None/0:03:00   0.755238
            api-o309vfqp3i58303hz3y1dk3g5   c2-standard-4    False       standard   0/1            10                  True     24/07, 14:36:25      None/0:03:00   0.755238
    """)

    subparser.set_defaults(func=list_machine_groups)
