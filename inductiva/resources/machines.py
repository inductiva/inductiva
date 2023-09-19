"""MachineGroup class to manage Google Cloud resources."""
import time

from absl import logging

import inductiva
import inductiva.client.models
from inductiva import api
from inductiva.client.apis.tags import instance_api


class MachineGroup():
    """Class to manage Google Cloud resources."""

    def __init__(
        self,
        machine_type: str,
        num_machines: int = 1,
        spot: bool = False,
        disk_size_gb: int = 40,
        zone: str = "europe-west1-b",
        register: bool = True,
    ) -> None:
        """Create a MachineGroup object.

        Args:
            machine_type: The type of GC machine to launch. Ex: "e2-medium".
            Check https://cloud.google.com/compute/docs/machine-resource for
            more information about machine types.
            num_machines: The number of virtual machines to launch.
            spot: Whether to use spot machines.
            disk_size_gb: The size of the disk in GB, recommended min. is 40 GB.
            zone: The zone where the machines will be launched.
        """
        self.id = None
        self.name = None
        self.create_time = None
        self.machine_type = machine_type
        self.num_machines = num_machines
        self.spot = spot
        self.disk_size_gb = disk_size_gb
        self.zone = zone
        self._started = False

        # Set the API configuration that carries the information from the client
        # to the backend.
        self._api = instance_api.InstanceApi(api.get_client())
        self._estimated_cost = None

        if register:
            self._register_machine_group()

    def _register_machine_group(self):
        instance_group_config = inductiva.client.models.InstanceGroupCreate(
            machine_type=self.machine_type,
            num_instances=self.num_machines,
            spot=self.spot,
            disk_size_gb=self.disk_size_gb,
            zone=self.zone,
        )
        logging.info("Registering machine group configuration with the API.")
        resp = self._api.register_instance_group(body=instance_group_config)
        self.id = resp.body["id"]
        self.name = resp.body["name"]
        logging.info("Registered machine group:")
        logging.info(" > ID: %s", self.id)
        logging.info(" > Name: %s", self.name)

    @classmethod
    def from_api_response(cls, resp: dict):
        """Creates a MachineGroup object from an API response."""

        machine_group = cls(
            machine_type=resp["machine_type"],
            num_machines=int(resp["num_instances"]),
            spot=bool(resp["spot"]),
            disk_size_gb=resp["disk_size_gb"],
            zone=resp["zone"],
            register=False,
        )
        machine_group.id = resp["id"]
        machine_group.name = resp["name"]
        machine_group.create_time = resp["create_time"]
        machine_group._started = True

        return machine_group

    def start(self):
        """Starts a machine group."""
        if self._started:
            logging.info("Attempting to start a machine group already started.")
            return

        if self.id is None or self.name is None:
            logging.info("Attempting to start an unregistered machine group. "
                         "Make sure you have called the constructor without "
                         "`register=False`.")
            return

        request_body = \
            inductiva.client.models.InstanceGroup(
                id=self.id,
                name=self.name,
                machine_type=self.machine_type,
                num_instances=self.num_machines,
                spot=self.spot,
                disk_size_gb=self.disk_size_gb,
                zone=self.zone,
            )
        try:
            logging.info("Starting a machine group. "
                         "This may take up to a few minutes.")
            logging.info("Note that stopping this local process will not "
                         "interrupt the creation of the machine group. "
                         "Please wait...")
            self._started = True
            start_time = time.time()
            self._api.start_instance_group(body=request_body)
            creation_time_mins = (time.time() - start_time) / 60

            logging.info("Machine group successfully started in %.2f mins.",
                         creation_time_mins)
            self._log_machine_group_info()

        except inductiva.client.ApiException as api_exception:
            raise api_exception

    def terminate(self):
        """Terminates a machine group."""
        if not self._started or self.id is None or self.name is None:
            logging.info("Attempting to terminate an unstarted machine group.")
            return

        try:
            logging.info("Terminating machine group. "
                         "This may take up to a few minutes.")
            start_time = time.time()

            request_body = \
                inductiva.client.models.InstanceGroup(
                    id=str(self.id),
                    name=self.name,
                    machine_type=self.machine_type,
                    num_instances=self.num_machines,
                    spot=self.spot,
                    disk_size_gb=self.disk_size_gb,
                    zone=self.zone,
                )

            self._api.delete_instance_group(body=request_body)
            termination_time_mins = (time.time() - start_time) / 60
            logging.info(
                "Machine group of %s machines successfully "
                "terminated in %.2f mins.", self.num_machines,
                termination_time_mins)

        except inductiva.client.ApiException as api_exception:
            raise api_exception

    def _get_estimated_cost(self) -> float:
        if self._estimated_cost is not None:
            return self._estimated_cost
        instance_price = self._api.get_instance_price({
            "machine_type": self.machine_type,
            "zone": self.zone,
        })
        if self.spot:
            estimated_cost = instance_price.body[
                "preemptible_price"] * self.num_machines
        else:
            estimated_cost = instance_price.body[
                "on_demand_price"] * self.num_machines

        self._estimated_cost = float(round(estimated_cost, 3))
        return self._estimated_cost

    def estimate_cloud_cost(self) -> float:
        """Returns the estimated cost per hour of a machine group.

        Note that this is an estimate of the cost incurred in cloud resources,
        and is not binding of the actual price of the machine group.
        """
        #TODO: Contemplate disk size in the price.
        estimated_cost = self._get_estimated_cost()
        logging.info("Estimated cloud cost per hour: %s $/h", estimated_cost)
        return estimated_cost

    def status(self):
        """Returns the status of a machine group if it exists.

        Otherwise returns None"""
        if self.name is None:
            logging.info(
                "Attempting to get the status of an unregistered machine "
                "group.")
            return

        response = self._api.get_group_status({"name": self.name})

        if response.body == "notFound":
            logging.info("Machine group does not exist: %s.", self.name)
        return response.body

    def _log_machine_group_info(self):
        """Logs the machine group info."""

        logging.info("Name: %s", self.name)
        logging.info("Machine type: %s", self.machine_type)
        logging.info("Number of machines: %s", self.num_machines)
        logging.info("Spot: %s", self.spot)
        logging.info("Disk size: %s GB", self.disk_size_gb)
        logging.info("Estimated cost per hour: %s $/h",
                     self._get_estimated_cost())
