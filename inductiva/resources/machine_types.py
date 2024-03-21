"""Available machine types and their number of cores."""

import json

import inductiva
from inductiva.client import ApiException
from inductiva.client.apis.tags import compute_api

from inductiva.resources import machines_base


def list_available_machines(provider: str):
    """List all available machines types."""

    resources_available = get_available_machine_types(provider)
    machine_types = []

    for machine in enumerate(resources_available):
        machine_types.append(machine["machine_type"])

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
