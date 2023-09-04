"""MachineGroup class to manage Google Cloud resources."""
import typing
import time

from absl import logging

import inductiva
from inductiva import api
from inductiva.client.apis.tags import instance_api


class MachineGroup():
    """Class to manage Google Cloud resources."""

    def __init__(
        self,
        machine_type: str,
        num_machines: int = 1,
        spot: bool = False,
        disk_size_gb: int = 30,
        zone: typing.Optional[str] = "europe-west1-b",
    ) -> None:
        """Create a MachineGroup object.

        Args:
            machine_type: The type of GC machine to launch. Ex: "e2-medium".
            Check https://cloud.google.com/compute/docs/machine-resource for
            more information about machine types.
            num_machines: The number of virtual machines to launch.
            spot: Whether to use spot machines.
            disk_size_gb: The size of the disk in GB, recommended min. is 30 GB.
            zone: The zone where the machines will be launched.
        """
        self.id = None
        self.name = None
        #TODO: Check if machine type is valid.
        self.machine_type = machine_type
        self.num_machines = num_machines
        self.spot = spot
        self.disk_size_gb = disk_size_gb
        self.zone = zone

        # Set the API configuration that carries the information from the client
        # to the backend.
        self._api = instance_api.InstanceApi(api.get_client())

    def start(self):
        """Starts a machine group."""
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
                         "This may take up to a few minutes.")
            start_time = time.time()
            instance_group = self._api.create_instance_group(
                body=instance_group_config)
            creation_time_mins = (time.time() - start_time) / 60

            self.id = instance_group.body["id"]
            self.name = instance_group.body["name"]
            self.estimated_price = self.estimate_price()

            logging.info("Machine group successfully created in %.2f mins.",
                         creation_time_mins)
            self._log_machine_group_info()

        except inductiva.client.ApiException as api_exception:
            raise api_exception

    def terminate(self):
        """Terminates a machine group."""

        try:
            logging.info("Terminating machine group."
                         "This may take up to a few minutes.")
            start_time = time.time()

            instance_group_config = \
                inductiva.client.model.instance_group.InstanceGroup(
                    name=self.name,
                    machine_type=self.machine_type,
                    num_instances=self.num_machines,
                    spot=self.spot,
                    disk_size_gb=self.disk_size_gb,
                    zone=self.zone)
            self._api.delete_instance_group(body=instance_group_config)
            termination_time_mins = (time.time() - start_time) / 60
            logging.info(
                "Machine group of %s machines successfully"
                "terminated in %.2f mins.", self.num_machines,
                termination_time_mins)

        except inductiva.client.ApiException as api_exception:
            raise api_exception

    def estimate_price(self):
        """Returns an estimated price per hour of a machine group."""
        #TODO: Contemplate disk size in the price.
        instance_price = self._api.get_instance_price({
            "machine_type": self.machine_type,
            "zone": self.zone,
        })
        if self.spot:
            estimated_price = instance_price.body[
                "preemptible_price"] * self.num_machines
        else:
            estimated_price = instance_price.body[
                "on_demand_price"] * self.num_machines
        estimated_price = float(round(estimated_price, 3))
        logging.info("Estimated price per hour: %s $/h", estimated_price)

        return estimated_price

    def status(self):
        """Returns the status of a machine group if it exists.

        Otherwise returns None"""
        response = self._api.get_group_status({"name": self.name})

        if response.body == "notFound":
            logging.info("Machine group does not exist: %s.", self.name)
        return response.body

    def _log_machine_group_info(self):
        """Logs the machine group info."""

        logging.info("ID: %s", self.id)
        logging.info("Machine type: %s", self.machine_type)
        logging.info("Number of machines: %s", self.num_machines)
        logging.info("Spot: %s", self.spot)
        logging.info("Disk size: %s GB", self.disk_size_gb)
        logging.info("Estimated cost per hour: %s $/h", self.estimated_price)
