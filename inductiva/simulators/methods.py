"""Methods to interact with the available simulators."""
from typing import Any
import inductiva
import json
from inductiva.client.apis.tags.simulators_api import SimulatorsApi


def list_available_images() -> Any:
    """Fetch the list of available simulator images from the API."""
    api = SimulatorsApi(inductiva.get_client())
    response = api.list_available_images(skip_deserialization=True).response
    return json.loads(response.data.decode("utf-8"))
