"""Methods related to managing executers and resource pools."""
from typing import Optional
from uuid import UUID

import inductiva
from inductiva import api
from inductiva.resources import Machines
from inductiva.client import ApiClient, ApiException
from inductiva.client.apis.tags.executers_api import ExecutersApi


def create_resource_pool() -> UUID:
    api_config = api.validate_api_key(inductiva.api_key)

    with ApiClient(api_config) as client:
        api_instance = ExecutersApi(client)

        try:
            api_response = api_instance.create_resource_pool()
        except ApiException as e:
            raise e

    return UUID(api_response.body["id"])


def launch_machine_group(
    machine_type: str,
    executer: str,
    num_machines: int = 1,
    spot: bool = True,
    disk_size_gb: int = 30,
    label: Optional[str] = None,
) -> None:

    resource_pool_id = create_resource_pool()

    return Machines(machine_type=machine_type,
                    num_machines=num_machines,
                    spot=spot,
                    executer=executer,
                    disk_size_gb=disk_size_gb,
                    label=label,
                    resource_pool_id=resource_pool_id)
