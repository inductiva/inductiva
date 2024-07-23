"""Available machine types and their number of cores."""
from dataclasses import dataclass
from typing import List, Union
import json

import inductiva
from inductiva.utils import format_utils
from inductiva.client import ApiException
from inductiva.client.apis.tags import compute_api


class ProviderType(format_utils.CaseInsensitiveEnum):
    """Enum to represent the provider of the machine to be launched."""
    GCP = "GCP"
    ICE = "ICE"
    LOCAL = "LOCAL"


@dataclass(frozen=True)
class MachineTypeInfo:
    """Dataclass to represent a machine type.
    This class has all the information related to a machine type.
    """

    machine_type: str
    provider_id: str
    num_cpus: int
    price: float
    ram_gb: int
    region: str
    spot: bool
    zone: str

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
        machine_family: str = None) -> List[MachineTypeInfo]:
    """Get all available machine types for a given provider.

    Args:
        provider (str): The provider to list the available machine types.
        machine_family (str): The machine family to filter the specific
            available machine types.
    Returns:
        list (MachineTypeInfo): List of available machine types.
    """

    provider = ProviderType(provider)

    api_client = compute_api.ComputeApi(inductiva.api.get_client())

    query_params = {"provider_id": provider.value}

    if machine_family is not None:
        query_params["machine_family"] = machine_family

    try:
        resp = api_client.list_available_machine_types(query_params).response

        response_body = json.loads(resp.data.decode("utf-8"))
        return [
            MachineTypeInfo(**machine_type) for machine_type in response_body
        ]
    except ApiException as e:
        raise e
