"""Methods related to managing executers and resource pools."""
from typing import Optional
from uuid import UUID

import inductiva
from inductiva import api
from inductiva.client import ApiClient, ApiException
from inductiva.client.apis.tags.executers_api import ExecutersApi
from inductiva.client.apis.tags.instance_api import InstanceApi
from inductiva.client.model.instance import Instance
from inductiva.client.model.instance_create import InstanceCreate
from inductiva.client.model.instance_group import InstanceGroup


def create_resource_pool() -> UUID:
    api_config = api.validate_api_key(inductiva.api_key)

    with ApiClient(api_config) as client:
        api_instance = ExecutersApi(client)

        try:
            api_response = api_instance.create_resource_pool()
        except ApiException as e:
            raise e

    return UUID(api_response.body["id"])


def launch_executer(
    name: str,
    machine_type: str,
    executer_type: str,
    spot: bool = True,
    resource_pool_id: Optional[UUID] = None,
    disk_size_gb: int = 10,
) -> None:

    api_config = api.validate_api_key(inductiva.api_key)

    with ApiClient(api_config) as client:
        api_instance = InstanceApi(client)

        body = InstanceCreate(
            name=name,
            machine_type=machine_type,
            image_name=executer_type,
            spot=spot,
            resource_pool_id=resource_pool_id,
            disk_size_gb=disk_size_gb,
        )

        try:
            api_instance.create_instance(body=body)
        except ApiException as e:
            raise e


def kill_executer(name: str) -> None:
    api_config = api.validate_api_key(inductiva.api_key)

    with ApiClient(api_config) as client:
        api_instance = InstanceApi(client)
        try:
            api_instance.delete_instance(body=Instance(name=name))
        except ApiException as e:
            raise e


def launch_executer_group(
    num_executers: int,
    name: str,
    machine_type: str,
    executer_type: str,
    spot: bool = True,
    resource_pool_id: Optional[UUID] = None,
    disk_size_gb: int = 10,
) -> None:

    api_config = api.validate_api_key(inductiva.api_key)

    with ApiClient(api_config) as client:
        api_instance = InstanceApi(client)

        body = InstanceGroup(
            name=name,
            machine_type=machine_type,
            image_name=executer_type,
            spot=spot,
            resource_pool_id=resource_pool_id,
            num_instances=num_executers,
            disk_size_gb=disk_size_gb,
        )

        try:
            api_instance.create_instance_group(body=body)
        except ApiException as e:
            raise e


def kill_executer_group(name: str) -> None:
    api_config = api.validate_api_key(inductiva.api_key)

    with ApiClient(api_config) as client:
        api_instance = InstanceApi(client)
        try:
            api_instance.delete_instance_group(body=Instance(name=name))
        except ApiException as e:
            raise e
