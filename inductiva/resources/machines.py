"""Containes Machines class to manage resources."""

from absl import logging
from typing import Optional
from uuid import UUID
import uuid
import inductiva
from inductiva import api

from inductiva.client import ApiClient, ApiException
from inductiva.client.apis.tags.instance_api import InstanceApi
from inductiva.client.model.instance import Instance
from inductiva.client.model.instance_group import InstanceGroup


class Machines():
    """Class to Google Cloud manage resources."""

    def __init__(
        self,
        machine_type: str,
        executer: str,
        num_machines: int = 1,
        spot: bool = False,
        disk_size_gb: int = 30,
        label: str = None,
        resource_pool_id: Optional[UUID] = None,
    ) -> None:
        """Initialize the Machines class.

        Args:
            machine_type: The type of GC machine to launch. Ex: "e2-medium".
            executer: The type of executer to launch. Ex: "gromacs", "openfoam".
            num_machines: The number of machines to launch.
            spot: Whether to use spot instances.
            disk_size_gb: The size of the disk in GB. (min. 30 GB)
            label: The label to assign to the machine group.
            resource_pool_id: The resource pool ID to use. A new resource pool
              id can be created using the
              inductiva.resources.create_resource_pool() method.
        """
        self.machine_type = machine_type
        self.num_machines = num_machines
        self.spot = spot
        self.disk_size_gb = disk_size_gb
        self.label = label
        self.executer = executer
        self.resource_pool_id = resource_pool_id
        self.name = self._generate_instance_name()

        self.api_config = api.validate_api_key(inductiva.api_key)

    def start(self):
        with ApiClient(self.api_config) as client:
            api_instance = InstanceApi(client)

            body = InstanceGroup(
                name=self.name,
                machine_type=self.machine_type,
                num_instances=self.num_machines,
                executer=self.executer,
                spot=self.spot,
                resource_pool_id=self.resource_pool_id,
                disk_size_gb=self.disk_size_gb,
            )
            try:
                logging.info("Creating a machine group. \
                             This may take a few minutes.")
                api_instance.create_instance_group(body=body)
                logging.info("Machine group is successfully created.")
            except ApiException as e:
                raise e

    def terminate(self):

        with ApiClient(self.api_config) as client:
            api_instance = InstanceApi(client)

            try:
                logging.info("Terminating machine group. \
                             This may take a few minutes.")
                api_instance.delete_instance_group(body=Instance(
                    name=self.name))
                logging.info("Machine group is successfully terminated.")

            except ApiException as e:
                raise e

    def estimate_cost(self):

        with ApiClient(self.api_config) as client:
            api_instance = InstanceApi(client)

            try:
                instance_price = api_instance.get_instance_price(body=Instance(
                    name=self.machine_type))
                logging.info("Estimated on-demand cost: %f",
                             instance_price.body["on_demand"])
                logging.info("Estimated spot cost: %f",
                             instance_price.body["preemptible"])
            except ApiException as e:
                raise e

    def _generate_instance_name(self):
        unique_id = uuid.uuid4().hex[:8]
        instance_name = f"api-{unique_id}"
        return instance_name
