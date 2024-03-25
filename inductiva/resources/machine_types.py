"""Available machine types and their number of cores."""

import json

import inductiva
from inductiva.client import ApiException
from inductiva.client.apis.tags import compute_api

from inductiva.resources import machines_base


def get_available():
    """Gets the available machine types.
    
    Currently, this fetches the available machine types from a yaml file
    into a dictionary that can now be used to parse the information where
    needed."""

    with open(MACHINE_TYPES_FILE, "r", encoding="utf-8") as file:
        machine_types = yaml.safe_load(file)

    return machine_types


def list_available_machines(provider: str):
    """List all available machines types."""

    resources_available = get_available()
    provider_resources = resources_available[provider]
    machine_types = []

    # Fetch the available CPU series for the given provider
    for cpu_series, series_info in provider_resources["cpu-series"].items():
        # Fetch the available RAM types and vCPUs info
        for ram_type, type_info in series_info["types"].items():
            vcpus = type_info["vcpus"]
            machine_types.extend(
                [f"{cpu_series}-{ram_type}-{vcpu}" for vcpu in vcpus])
            if type_info["lssd"]:
                for vcpu in vcpus:
                    machine_types.append(f"{cpu_series}-{ram_type}-{vcpu}-lssd")

    return tuple(machine_types)


def get_available_machine_types(provider: str, machine_family: str = None):
    """Get all available machine types for a given provider.
    
    Args:
        provider (str): The provider to list the available machine types.
        machine_family (str): The machine family to filter the specific
            available machine types."""

    provider = machines_base.ProviderType(provider.upper())

    api_client = compute_api.ComputeApi(inductiva.api.get_client())

    query_params = {"provider_id": provider.value}

    if machine_family is not None:
        query_params["machine_family"] = machine_family

    try:
        resp = api_client.list_available_machine_types(query_params).response

        response_body = json.loads(resp.data.decode("utf-8"))

        return response_body
    except ApiException as e:
        raise e
