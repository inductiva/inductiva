"""Containes MachineGroup class to manage GC resources."""
import typing
import uuid
import time

from absl import logging

import inductiva


class MachineGroup():
    """Class to manage Google Cloud resources."""

    def __init__(
        self,
        machine_type: str,
        num_machines: int = 1,
        spot: bool = False,
        disk_size_gb: int = 20,
        zone: typing.Optional[str] = "europe-west1-b",
    ) -> None:
        """Create a MachineGroup object.

        Args:
            machine_type: The type of GC machine to launch. Ex: "e2-medium".
            num_machines: The number of machines to launch.
            spot: Whether to use spot machines.
            disk_size_gb: The size of the disk in GB, recommended min. is 20 GB.
            name: The name to assign to the machine group.
        """
        #TODO: Check if machine type is valid.
        self.machine_type = machine_type
        self.num_machines = num_machines
        self.spot = spot
        self.disk_size_gb = disk_size_gb
        #TODO: Pass the name generation to the backend
        self.name = self._generate_instance_name()
        self.zone = zone
        self.estimated_cost = -1  # The negative implies that the cost is unknown.

        # Set the API configuration that carries the information from the client
        # to the backend.
        self.api_config = inductiva.api.validate_api_key(inductiva.api_key)

    def start(self):
        """Starts a machine group."""

        with inductiva.client.ApiClient(self.api_config) as client:
            api_instance = inductiva.client.apis.tags.instance_api.InstanceApi(
                client)

            #TODO: Set this creation on the backend
            self.machine_group_id = \
                inductiva.resources.utils.create_machine_group_id()
            instance_group_config = inductiva.client.model.instance_group.InstanceGroup(
                name=self.name,
                machine_type=self.machine_type,
                num_instances=self.num_machines,
                spot=self.spot,
                resource_pool_id=self.machine_group_id,
                disk_size_gb=self.disk_size_gb,
                zone=self.zone,
            )
            try:
                logging.info("Creating a machine group. \
                             This may take a few minutes.")
                start_time = time.time()
                #TODO: Receive here the machine name and ID.
                api_instance.create_instance_group(body=instance_group_config)
                creation_time_mins = (time.time() - start_time) / 60

                self.estimated_cost = self._compute_estimated_cost(api_instance)
                logging.info(
                    "Machine group with the specified settings successfully \
                        created.")
                self._log_machine_group_info()
                logging.info("Machine group successfully created in %s mins.",
                             creation_time_mins)

            except inductiva.client.ApiException as api_exception:
                raise api_exception

    def terminate(self):
        """Terminates a machine group."""

        with inductiva.client.ApiClient(self.api_config) as client:
            api_instance = inductiva.client.apis.tags.instance_api.InstanceApi(
                client)

            try:
                logging.info("Terminating machine group. \
                             This may take a few minutes.")
                start_time = time.time()
                api_instance.delete_instance_group(
                    body=inductiva.client.model.instance.Instance(
                        name=self.name))
                termination_time_mins = (time.time() - start_time) / 60
                logging.info(
                    "Machine group of %s machines successfully \
                    terminated in % s mins.", self.num_machines,
                    termination_time_mins)

            except inductiva.client.ApiException as api_exception:
                raise api_exception

    def _compute_estimated_cost(self, api_instance):
        """Returns an estimated cost per hour of a machine group."""

        instance_price = api_instance.get_instance_price(
            body=inductiva.client.model.instance.Instance(
                name=self.machine_type, zone=self.zone))

        if self.spot:
            estimated_price = instance_price.body[
                "on_demand"] * self.num_machines
        else:
            estimated_price = instance_price.body[
                "preemptible"] * self.num_machines

        return estimated_price

    def _generate_instance_name(self):
        """Generate instance name based on unique_id.
        
        TODO: Pass this to the backend."""

        unique_id = uuid.uuid4().hex[:8]
        instance_name = f"api-{unique_id}"
        return instance_name

    def _log_machine_group_info(self):
        """Logs the machine group info."""

        logging.info("Machine type: %s", self.machine_type)
        logging.info("Number of machines: %s", self.num_machines)
        logging.info("Spot: %s", self.spot)
        logging.info("Disk size: %s GB", self.disk_size_gb)
        logging.info("Machine group ID: %s", self.machine_group_id)
        logging.info("Estimated cost per hour: %s $/h", self.estimated_cost)
        logging.info("Name: %s", self.name)
