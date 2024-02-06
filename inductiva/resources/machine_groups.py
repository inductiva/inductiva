"""Functions to manage or retrieve user resources."""
from typing import Optional
from absl import logging
import inductiva
import inductiva.client
from inductiva.client.apis.tags import compute_api
from inductiva.utils import format_utils
from inductiva import resources
from . import machines_base


def estimate_machine_cost(machine_type: str, spot: bool = False):
    """Estimate the cloud cost of one machine per hour in US dollars.

    Args:
        machine_type: The type of GC machine to launch. Ex: "e2-standard-4".
            Check https://cloud.google.com/compute/docs/machine-resource for
            more information about machine types.
        zone: The zone where the machines will be launched.
        spot: Whether to use spot machines.
    """

    api = compute_api.ComputeApi(inductiva.api.get_client())

    instance_price = api.get_instance_price({
        "machine_type": machine_type,
    })

    if spot:
        estimated_cost = instance_price.body["preemptible_price"]
    else:
        estimated_cost = instance_price.body["on_demand_price"]

    return float(estimated_cost)


def _machine_group_list_to_str(machine_group_list) -> str:
    """Returns a string representation of a list of machine groups."""
    columns = [
        "Name",
        "Machine Type",
        "Elastic",
        "Type",
        "# machines",
        "Data Size in GB",
        "Spot",
        "Started at (UTC)",
    ]
    rows = []

    for machine_group in machine_group_list:
        is_elastic = False
        resource_type = machines_base.ResourceType.STANDARD.value
        spot = machine_group.spot if hasattr(machine_group, "spot") else False

        if isinstance(machine_group, resources.ElasticMachineGroup):
            num_active_machines = machine_group.num_active_machines
            max_machines = machine_group.max_machines
            num_active_machines = f"{num_active_machines}/{max_machines}"
            is_elastic = True
        else:
            num_active_machines = machine_group.num_machines
            if isinstance(machine_group, resources.MPICluster):
                resource_type = machines_base.ResourceType.MPI.value

        rows.append([
            machine_group.name, machine_group.machine_type, is_elastic,
            resource_type, num_active_machines, machine_group.data_disk_gb,
            spot, machine_group.create_time
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


def _fetch_machine_groups_from_api():
    """Get all active machine groups of a user from the API."""
    try:
        api = compute_api.ComputeApi(inductiva.api.get_client())
        response = api.list_active_user_instance_groups()
        if len(response.body) == 0:
            print("No active computational resources found.")
            return response.body

        return response.body

    except inductiva.client.ApiException as api_exception:
        raise api_exception


# pylint: disable=redefined-builtin
def list():
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

    machine_group_list = get()
    if len(machine_group_list) != 0:
        print("Active Resources:")
        print(_machine_group_list_to_str(machine_group_list))


def get():
    """Returns a list of 'Resource' objects."""

    # Retrive the active resource names
    machine_groups = _fetch_machine_groups_from_api()
    machine_group_list = []

    for mg in machine_groups:
        if mg["is_elastic"]:
            mg_class = resources.ElasticMachineGroup
        elif mg["type"] == "standard":
            mg_class = resources.MachineGroup
        elif mg["type"] == "mpi":
            mg_class = resources.MPICluster
        else:
            raise ValueError("Unknown resource configuration.")
        machine_group_list.append(mg_class.from_api_response(mg))

    return machine_group_list


def get_cheapest_machine_type(num_cpus: int,
                              ram_gb: Optional[int] = None,
                              spot: bool = False):
    """Get the machine type with the lowest price.

    If the machine type with the exact requirements is not found, the closest
    machine type with higher resources will be returned. For example, if a
    machine type with 5 CPUs and 20 GB of RAM is requested, but the closest
    machine type has 8 CPUs and 32 GB of RAM, the latter will be returned.

    Args:
        num_cpus: The minimum number of CPUs required.
        ram_gb: The minimum RAM in GB required. If None, the minimum RAM for
          the given number of CPUs will be used.
        spot: Whether the machine is spot or not.
    """
    spot = "t" if spot else "f"
    api = compute_api.ComputeApi(inductiva.api.get_client())
    body = {"num_cpus": num_cpus, "spot": spot}
    if ram_gb:
        body["ram_gb"] = ram_gb
    response = api.get_machine_type(body)

    logging.info("Machine type: %s", response.body["machine_type"])
    logging.info("CPUs: %s", response.body["num_cpus"])
    logging.info("RAM in GB: %s", response.body["ram_gb"])
    logging.info("Estimated cost per hour: %s $/h",
                 round(response.body["price"], 5))

    return response.body["machine_type"]
