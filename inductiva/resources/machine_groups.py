"""Functions to manage or retrieve user resources."""
from collections import defaultdict

import inductiva
import inductiva.client
from inductiva.client.apis.tags import compute_api
from inductiva import resources


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


def get_machine_dict(machines):
    """Get a dictionary with the information of the machines."""
    column_names = [
        "Host Name",
        "Started",
        "Status",
        "Last seen",
        "Running task",
    ]
    table = defaultdict(list, {key: [] for key in column_names})
    for machine in machines:
        status = "Off" if isinstance(machine["terminated_at"],
                                     str) else "Active"
        task_id = machine["current_task_id"] if isinstance(
            machine["current_task_id"], str) else None
        table["Host Name"].append(machine["host_name"])
        table["Started"].append(machine["started_at"])
        table["Status"].append(status)
        table["Last seen"].append(machine["last_seen_at"])
        table["Running task"].append(task_id)

    return table


def _fetch_machine_groups_from_api():
    """Get all active machine groups of a user from the API."""
    try:
        api = compute_api.ComputeApi(inductiva.api.get_client())
        response = api.list_active_user_instance_groups()

        return response.body

    except inductiva.client.ApiException as api_exception:
        raise api_exception


def get_by_name(machine_name: str):
    """Returns the machine group corresponding to `machine_name`."""
    try:
        api = compute_api.ComputeApi(inductiva.api.get_client())
        response = api.get_vm_group_by_name({"name": machine_name}).body
        mg_class = _get_machine_group_class(response["type"],
                                            response["is_elastic"])
        return mg_class.from_api_response(response)
    except inductiva.client.ApiException as api_exception:
        raise api_exception


def _get_machine_group_class(machine_type: str, is_elastic: bool):
    """Returns the class of the machine group"""
    if is_elastic:
        mg_class = resources.ElasticMachineGroup
    elif machine_type == "standard":
        mg_class = resources.MachineGroup
    elif machine_type == "mpi":
        mg_class = resources.MPICluster
    else:
        raise ValueError("Unknown resource configuration.")
    return mg_class


def get():
    """Returns a list of 'Resource' objects."""

    # Retrive the active resource names
    machine_groups = _fetch_machine_groups_from_api()
    machine_group_list = []

    for mg in machine_groups:
        mg_class = _get_machine_group_class(mg["type"], mg["is_elastic"])
        machine_group_list.append(mg_class.from_api_response(mg))

    return machine_group_list
