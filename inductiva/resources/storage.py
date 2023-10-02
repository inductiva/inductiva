"""Methods to interact with the user storage resources."""
import inductiva
from inductiva.client.apis.tags import instance_api


def get_space_used():
    """Returns the occupied storage size in GB."""
    try:
        api = instance_api.InstanceApi(inductiva.api.get_client())
        storage_size = api.get_storage_size()
        storage_size = float(round(storage_size.body, 3))
        print(f"Total storage used: {storage_size} GB")
        return storage_size
    except inductiva.client.ApiException as api_exception:
        raise api_exception

