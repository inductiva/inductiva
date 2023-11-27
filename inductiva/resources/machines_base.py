"""Base class for machine groups."""
import time

from absl import logging

import inductiva
import inductiva.client.models
from inductiva import api
from inductiva.client.apis.tags import compute_api


class BaseMachineGroup():
    """Base class to manage Google Cloud resources."""

    def __init__(self,
                 machine_type: str,
                 spot: bool = False,
                 disk_size_gb: int = 70,
                 zone: str = "europe-west1-b",
                 register: bool = True) -> None:
        """Create a BaseMachineGroup object.

        Args:
            machine_type: The type of GC machine to launch. Ex: "e2-standard-4".
              Check https://cloud.google.com/compute/docs/machine-resource for
              more information about machine types.
            spot: Whether to use spot machines.
            disk_size_gb: The size of the disk in GB, recommended min. is 60 GB.
            zone: The zone where the machines will be launched.
        """
        self.machine_type = machine_type
        self.spot = spot
        self.disk_size_gb = disk_size_gb
        self.zone = zone
        self._id = None
        self._name = None
        self.create_time = None
        self._started = False
        self.register = register

        # Set the API configuration that carries the information from the client
        # to the backend.
        self._api = compute_api.ComputeApi(api.get_client())
        self._estimated_cost = None

    @property
    def id(self):
        return self._id

    @property
    def name(self):
        return self._name

    def _register_machine_group(self, **kwargs):
        instance_group_config = inductiva.client.models.InstanceGroupCreate(
            machine_type=self.machine_type,
            spot=self.spot,
            disk_size_gb=self.disk_size_gb,
            zone=self.zone,
            **kwargs,
        )
        logging.info("Registering machine group configurations:")
        resp = self._api.register_instance_group(body=instance_group_config)
        self._id = resp.body["id"]
        self._name = resp.body["name"]
        self.register = False

    @classmethod
    def from_api_response(cls, resp: dict):
        """Creates a MachineGroup object from an API response."""

        machine_group = cls(
            machine_type=resp["machine_type"],
            spot=bool(resp["spot"]),
            disk_size_gb=resp["disk_size_gb"],
            zone=resp["zone"],
            register=False,
        )
        machine_group._id = resp["id"]
        machine_group._name = resp["name"]
        machine_group.create_time = resp["create_time"]
        machine_group._started = True

        return machine_group

    def start(self, **kwargs):
        """Starts a machine group.

        Args:
            **kwargs: Depending on the type of machine group to be started,
              this can be num_machines, max_machines, min_machines,
              and is_elastic."""
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
                spot=self.spot,
                disk_size_gb=self.disk_size_gb,
                zone=self.zone,
                **kwargs,
            )
        try:
            logging.info("Starting machine group. "
                         "This may take a few minutes.")
            logging.info("Note that stopping this local process will not "
                         "interrupt the creation of the machine group. "
                         "Please wait...")
            start_time = time.time()
            if kwargs.get("is_elastic", True):
                self._api.start_elastic_instance_group(body=request_body)
            else:
                self._api.start_instance_group(body=request_body)
            creation_time_mins = (time.time() - start_time) / 60
            self._started = True

            logging.info("Machine group successfully started in %.2f mins.\n",
                         creation_time_mins)

        except inductiva.client.ApiException as api_exception:
            raise api_exception

    def terminate(self, **kwargs):
        """Terminates a machine group."""
        if not self._started or self.id is None or self.name is None:
            logging.info("Attempting to terminate an unstarted machine group.")
            return

        try:
            logging.info("Terminating machine group. "
                         "This may take a few minutes.")
            start_time = time.time()

            request_body = \
                inductiva.client.models.InstanceGroup(
                    id=self.id,
                    name=self.name,
                    machine_type=self.machine_type,
                    spot=self.spot,
                    disk_size_gb=self.disk_size_gb,
                    zone=self.zone,
                    **kwargs,
                )

            self._api.delete_instance_group(body=request_body)
            termination_time_mins = (time.time() - start_time) / 60
            logging.info(
                "Machine group '%s' successfully "
                "terminated in %.2f mins.\n", self.name, termination_time_mins)

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
            estimated_cost = instance_price.body["preemptible_price"]
        else:
            estimated_cost = instance_price.body["on_demand_price"]

        self._estimated_cost = estimated_cost
        return self._estimated_cost

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

        logging.info("> Name: %s", self.name)
        logging.info("> Machine type: %s", self.machine_type)
        # TODO: Not yet available to users
        logging.info("> Spot: %s", self.spot)
        logging.info("> Disk size: %s GB", self.disk_size_gb)
