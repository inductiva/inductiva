"""Methods related to managing executers and resource pools."""
from uuid import UUID

import inductiva
from inductiva import api
from inductiva.client import ApiClient, ApiException
from inductiva.client.apis.tags.executers_api import ExecutersApi


def create_machine_group_id() -> UUID:
    api_config = api.validate_api_key(inductiva.api_key)

    with ApiClient(api_config) as client:
        api_instance = ExecutersApi(client)

        try:
            api_response = api_instance.create_resource_pool()
        except ApiException as e:
            raise e

    return UUID(api_response.body["id"])
