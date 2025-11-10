"""Methods to interact with the available simulators."""
from typing import Dict, List
import inductiva
import inductiva.machines_catalogue_client


def list_available_images() -> Dict[str, Dict[str, List[str]]]:
    """Fetch the list of available simulator images from the API."""
    client = inductiva.api.get_machines_catalogue_client()
    api = inductiva.machines_catalogue_client.SimulatorsApi(client)
    response = api.list_available_images_simulators_available_images_get()
    return response.model_dump()
