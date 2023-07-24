"""Containes MachineGroup class to manage GC resources."""
from absl import logging
from typing import Optional
import uuid
import inductiva
from inductiva import api
from inductiva import resources

from inductiva.client import ApiClient, ApiException
from inductiva.client.apis.tags.instance_api import InstanceApi
from inductiva.client.model.instance import Instance
from inductiva.client.model.instance_group import InstanceGroup


class MachineGroup():
    """Class to manage Google Cloud resources."""

    def __init__(
        self,
        machine_type: str,
        num_machines: int = 1,
        spot: bool = False,
        disk_size_gb: int = 30,
        label: Optional[str] = None,
    ) -> None:
        """Create a MachineGroup object.

        Args:
            machine_type: The type of GC machine to launch. Ex: "e2-medium".
            num_machines: The number of machines to launch.
            spot: Whether to use spot instances.
            disk_size_gb: The size of the disk in GB, recommended min. is 30 GB.
            label: The label to assign to the machine group.
        """
        self.machine_type = machine_type
        self.num_machines = num_machines
        self.spot = spot
        self.disk_size_gb = disk_size_gb
        self.label = label
        self.machine_group_id = resources.create_machine_group_id()
        self.name = self._generate_instance_name()

        self.api_config = api.validate_api_key(inductiva.api_key)

    def start(self):
        with ApiClient(self.api_config) as client:
            api_instance = InstanceApi(client)

            body = InstanceGroup(
                name=self.name,
                machine_type=self.machine_type,
                num_instances=self.num_machines,
                spot=self.spot,
                resource_pool_id=self.machine_group_id,
                disk_size_gb=self.disk_size_gb,
            )
            try:
                logging.info("Creating a machine group. \
                             This may take a few minutes.")
                api_instance.create_instance_group(body=body)
                logging.info("Machine group successfully created.")
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
                logging.info("Machine group successfully terminated.")

            except ApiException as e:
                raise e

    def estimate_cost(self):
        """Returns an estimated cost per hour of a machine group."""
        with ApiClient(self.api_config) as client:
            api_instance = InstanceApi(client)

            try:
                instance_price = api_instance.get_instance_price(body=Instance(
                    name=self.machine_type))
                logging.info(
                    "Estimated cost of a machine group is %s$ per hour.",
                    self._get_cost(instance_price))
            except ApiException as e:
                raise e
        return self._get_cost(instance_price)

    def _get_cost(self, instance_price):
        if self.spot:
            return instance_price.body["on_demand"] * self.num_machines
        else:
            return instance_price.body["preemptible"] * self.num_machines

    def _generate_instance_name(self):
        unique_id = uuid.uuid4().hex[:8]
        instance_name = f"api-{unique_id}"
        return instance_name
