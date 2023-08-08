"""Methods related to managing executers and resource pools."""
import uuid

import inductiva


def create_machine_group_id() -> uuid.UUID:
    """Request the creation of a machine group ID."""

    api_config = inductiva.api.validate_api_key(inductiva.api_key)

    with inductiva.client.ApiClient(api_config) as client:
        api_instance = inductiva.client.apis.tags.executers_api.ExecutersApi(
            client)

        try:
            api_response = api_instance.create_resource_pool()
        except inductiva.client.ApiException as exception:
            raise exception

    return uuid.UUID(api_response.body["id"])
