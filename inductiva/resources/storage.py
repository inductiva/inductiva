"""Methods to interact with the user storage resources."""
import inductiva
from inductiva.client.apis.tags import instance_api


def get_space_used():
    """Returns the occupied storage size in GB."""
    try:
        api = instance_api.InstanceApi(inductiva.api.get_client())
        instance_price = api.get_storage_size()
        instance_price = float(round(instance_price.body, 3))
        print(f"Total storage used: {instance_price} GB")
        return instance_price
    except inductiva.client.ApiException as api_exception:
        raise api_exception
