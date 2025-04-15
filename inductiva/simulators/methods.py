"""Methods to interact with the available simulators."""
from typing import Any
import inductiva
from inductiva.client.apis.tags.simulators_api import SimulatorsApi






def list_available_images() -> Any:
    """Fetch the list of available simulator images from the API."""
    api = SimulatorsApi(inductiva.get_client())
    response =api.list_available_images()
    return response.body
