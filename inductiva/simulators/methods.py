"""Methods to interact with the available simulators."""
from typing import Dict, List
import inductiva
import inductiva.client
import json


def list_available_images() -> Dict[str, Dict[str, List[str]]]:
    """Fetch the list of available simulator images from the API."""
    api = inductiva.client.SimulatorsApi(inductiva.get_client())
    response = api.list_available_images_without_preload_content()
    return json.loads(response.data.decode("utf-8"))
