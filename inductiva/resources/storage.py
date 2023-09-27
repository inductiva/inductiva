"""Methods to interact with the user storage resources."""
import inductiva
from inductiva.client.apis.tags import instance_api


def get_space_used():
    """Returns the occupied storage size in GB."""
    try:
        api = instance_api.InstanceApi(inductiva.api.get_client())
        response = api.get_storage_size()
        sstorage_used = float(round(response.body, 3))
        print(f"Total storage used: {sstorage_used} GB")
        return sstorage_used
    except inductiva.client.ApiException as api_exception:
        raise api_exception
