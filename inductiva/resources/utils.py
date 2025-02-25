"""Functions to manage or retrieve user resources."""
from collections import defaultdict
from dataclasses import dataclass
import json
from typing import List, Optional, Tuple, Union

import inductiva
import inductiva.client
from inductiva.client.apis.tags import compute_api
from inductiva.client import ApiException
from inductiva.utils import format_utils


class ProviderType(format_utils.CaseInsensitiveEnum):
    """Enum to represent the provider of the machine to be launched."""
    GCP = "GCP"
    LOCAL = "LOCAL"


@dataclass(frozen=True)
class MachineTypeInfo:
    """Dataclass to represent a machine type.
    This class has all the information related to a machine type.
    """

    threads_per_core: int
    machine_type: str
    provider_id: str
    num_cpus: int
    price: float
    ram_gb: int
    region: str
    spot: bool
    zone: str
    num_gpus: Optional[int]
    gpu_name: Optional[str]

    def __getitem__(self, item):
        return getattr(self, item)


def list_available_machines(provider: Union[str, ProviderType]):
    """List all available machines types."""

    resources_available = get_available_machine_types(provider)
    machine_types = []

    for machine in resources_available:
        machine_types.append(machine["machine_type"])

    return tuple(machine_types)


def get_available_machine_types(
    provider: Union[str, ProviderType] = ProviderType.GCP,
    machine_families: Optional[List[str]] = None,
    machine_configs: Optional[List[str]] = None,
    vcpus_range: Optional[Tuple[int, int]] = None,
    memory_range: Optional[Tuple[int, int]] = None,
    price_range: Optional[Tuple[float, float]] = None,
    spot: Optional[bool] = None,
) -> List[MachineTypeInfo]:
    """Get all available machine types from a specified provider,
    allowing for multiple filtering parameters.
    
    Args:
        provider (Union[str, ProviderType]): 
            The cloud provider for which to list available machine types. 
            Defaults to `ProviderType.GCP`.
        machine_families (Optional[List[str]]): 
            A list of machine families to filter the available machine types
            (e.g., "c2", "c2d", "n4", etc).
        machine_configs (Optional[List[str]]): 
            A list of specific machine configurations to filter the results
            (e.g., "highcpu", "highmem", "standard").
        vcpus_range (Optional[Tuple[int, int]]): 
            A tuple defining the range of virtual CPUs (vCPUs) to filter the 
            machine types.
        memory_range (Optional[Tuple[int, int]]): 
            A tuple defining the range of memory (in MB) to filter the machine 
            types.
        price_range (Optional[Tuple[float, float]]): 
            A tuple defining the price range (in the respective currency) to 
            filter the machine types.
        spot (Optional[bool]): 
            If set to True, filters for spot instances; if False, filters for 
            on-demand instances.

    Returns:
        List[MachineTypeInfo]:
            A list of available machine types that match the specified criteria.
    """

    provider = ProviderType(provider)

    api_client = compute_api.ComputeApi(inductiva.api.get_client())

    query_params = {"provider_id": provider.value}

    if machine_families:
        query_params["machine_families"] = machine_families

    if machine_configs:
        query_params["machine_configs"] = machine_configs

    if vcpus_range:
        query_params["vcpus_range"] = vcpus_range

    if memory_range:
        query_params["memory_range"] = memory_range

    if price_range:
        query_params["price_range"] = price_range

    if spot is not None:
        query_params["spot"] = str(spot).lower()

    try:
        resp = api_client.list_available_machine_types(query_params).response

        response_body = json.loads(resp.data.decode("utf-8"))
        return [
            MachineTypeInfo(**machine_type) for machine_type in response_body
        ]
    except ApiException as e:
        raise e


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
