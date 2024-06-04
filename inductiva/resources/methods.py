"""Methods related with all resource types."""
from typing import Any, Dict, List

from inductiva.resources.machine_types import get_available_machine_types


def query(provider: str = "GCP", query_filter=None) -> List[Dict[str, Any]]:
    """Query available resources.

    This function returns a list of available resources for a given provider. And
    optionally filters the results based on a query_filter. The query_filter
    should be a callable that takes a single argument, the machine type (str).

    Args:
        provider (str): Provider name.
        query_filter (callable): Filter function.
            The resulting list will contain only the elements for which the
            callable returns True.
            The callable should take a single argument, the machine type (str).
    returns:
        List[Dict[str, Any]]: List of available resources.
    """
    queried_machines = []

    all_machines = get_available_machine_types(provider)

    if query_filter is None:
        return all_machines
    for machine in all_machines:
        if callable(query_filter):
            if query_filter(machine["machine_type"]):
                queried_machines.append(machine)
        else:
            raise ValueError("query_filter must be a callable")
    return queried_machines
