"""MachineGroup class to manage Google Cloud resources."""
import typing
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
            Check https://cloud.google.com/compute/docs/machine-resource for 
            more information about machine types.
            num_machines: The number of virtual machines to launch.
            spot: Whether to use spot machines.
            disk_size_gb: The size of the disk in GB, recommended min. is 20 GB.
            zone: The zone where the machines will be launched.
        """
        self.id = None
        #TODO: Check if machine type is valid.
        self.machine_type = machine_type
        self.num_machines = num_machines
        self.spot = spot
        self.disk_size_gb = disk_size_gb
        self.zone = zone
        # The negative imply that the cost is unknown.
        self.estimated_price = -1

        # Set the API configuration that carries the information from the client
        # to the backend.
        self.api_config = inductiva.api.validate_api_key(inductiva.api_key)

    def start(self):
        """Starts a machine group."""

        with inductiva.client.ApiClient(self.api_config) as client:
            api_instance = inductiva.client.apis.tags.instance_api.InstanceApi(
                client)

            instance_group_config = \
                inductiva.client.model.instance_group.InstanceGroup(
                    machine_type=self.machine_type,
                    num_instances=self.num_machines,
                    spot=self.spot,
                    disk_size_gb=self.disk_size_gb,
                    zone=self.zone,
                )
            try:
                logging.info("Creating a machine group."
                             "This may take a few minutes.")
                start_time = time.time()
                instance_group = api_instance.create_instance_group(
                    body=instance_group_config)
                creation_time_mins = time.time() - start_time

                self.id = instance_group.body["id"]
                self.name = instance_group.body["name"]
                #self.estimated_price = self._compute_estimated_price(
                #    api_instance)

                logging.info("Machine group successfully created in %.2f s.",
                             creation_time_mins)
                self._log_machine_group_info()

            except inductiva.client.ApiException as api_exception:
                raise api_exception

    def terminate(self):
        """Terminates a machine group."""

        with inductiva.client.ApiClient(self.api_config) as client:
            api_instance = inductiva.client.apis.tags.instance_api.InstanceApi(
                client)

            instance_group_config = \
                inductiva.client.model.instance_group.InstanceGroup(
                    name=self.name,
                    machine_type=self.machine_type,
                    num_instances=self.num_machines,
                    spot=self.spot,
                    disk_size_gb=self.disk_size_gb,
                    zone=self.zone,
                )
            try:
                logging.info("Terminating machine group."
                             "This may take a few minutes.")
                start_time = time.time()
                api_instance.delete_instance_group(body=instance_group_config)
                termination_time_mins = time.time() - start_time
                logging.info(
                    "Machine group of %s machines successfully "
                    "terminated in %.2f s.", self.num_machines,
                    termination_time_mins)

            except inductiva.client.ApiException as api_exception:
                raise api_exception

    def _compute_estimated_price(self, api_instance):
        """Returns an estimated price per hour of a machine group."""
        #TODO: Contemplate disk size in the price.
        instance_price = api_instance.get_instance_price(
            body=inductiva.client.model.instance.Instance(id=self.machine_type,
                                                          zone=self.zone))

        if self.spot:
            estimated_price = instance_price.body[
                "on_demand"] * self.num_machines
        else:
            estimated_price = instance_price.body[
                "preemptible"] * self.num_machines

        return estimated_price

    def _log_machine_group_info(self):
        """Logs the machine group info."""

        logging.info("ID: %s", self.id)
        logging.info("Machine type: %s", self.machine_type)
        logging.info("Number of machines: %s", self.num_machines)
        logging.info("Spot: %s", self.spot)
        logging.info("Disk size: %s GB", self.disk_size_gb)
        # We won't offer price values just yet, due to a bug.
        # TODO: Fix the price computation.
        # logging.info("Estimated cost per hour: %s $/h", self.estimated_price)
