"""Base class for machine groups."""
import time
import enum

from absl import logging

import inductiva
import inductiva.client.models
from inductiva import api
from inductiva.client.apis.tags import compute_api


class ResourceType(enum.Enum):
    """Enum to represent the type of machine to be launched."""

    STANDARD = "standard"
    MPI = "mpi"


class BaseMachineGroup():
    """Base class to manage Google Cloud resources."""

    def __init__(self,
                 machine_type: str,
                 disk_size_gb: int = 70,
                 register: bool = True) -> None:
        """Create a BaseMachineGroup object.

        Args:
            machine_type: The type of GC machine to launch. Ex: "e2-standard-4".
              Check https://cloud.google.com/compute/docs/machine-resource for
              more information about machine types.
            spot: Whether to use spot machines.
            disk_size_gb: The size of the disk in GB, recommended min. is 60 GB.
            register: Bool that indicates if a machine group should be register
                or if it was already registered. If set to False by users on
                initialization, then, the machine group will not be able to be
                started. This serves has an helper argument for retrieving
                already registered machine groups that can be started, for
                example, when retrieving with the `machines_groups.get` method.
                Users should not set this argument in anyway.
        """
        allowed_machine_types = [
            machine_type.name + str(core)
            for machine_type in inductiva.resources.machine_types.MachineType
            for core in machine_type.value
        ]
        assert machine_type in allowed_machine_types,\
            "Machine type not supported."
        self.machine_type = machine_type
        self.disk_size_gb = disk_size_gb
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
        """Register machine group configuration in API.
        
        Returns:
            The unique ID and name identifying the machine on the API."""

        instance_group_config = inductiva.client.models.GCPVMGroup(
            machine_type=self.machine_type,
            disk_size_gb=self.disk_size_gb,
            **kwargs,
        )
        logging.info("Registering machine group configurations:")
        resp = self._api.register_vm_group(body=instance_group_config)
        self._id = resp.body["id"]
        self._name = resp.body["name"]
        self.register = False
        self._log_machine_group_info()

    @classmethod
    def from_api_response(cls, resp: dict):
        """Creates a MachineGroup object from an API response."""

        machine_group = cls(
            machine_type=resp["machine_type"],
            disk_size_gb=resp["disk_size_gb"],
            register=False,
        )
        machine_group._id = resp["id"]
        machine_group._name = resp["name"]
        machine_group.create_time = resp["creation_timestamp"]
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
            inductiva.client.models.GCPVMGroup(
                id=self.id,
                name=self.name,
                machine_type=self.machine_type,
                disk_size_gb=self.disk_size_gb,
                **kwargs,
            )
        try:
            logging.info("Starting machine group. "
                         "This may take a few minutes.")
            logging.info("Note that stopping this local process will not "
                         "interrupt the creation of the machine group. "
                         "Please wait...")
            start_time = time.time()
            self._api.start_vm_group(body=request_body)
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
                inductiva.client.models.GCPVMGroup(
                    id=self.id,
                    name=self.name,
                    machine_type=self.machine_type,
                    disk_size_gb=self.disk_size_gb,
                    **kwargs,
                )

            self._api.delete_vm_group(body=request_body)
            termination_time_mins = (time.time() - start_time) / 60
            logging.info(
                "Machine group '%s' successfully "
                "terminated in %.2f mins.\n", self.name, termination_time_mins)

        except inductiva.client.ApiException as api_exception:
            raise api_exception

    def _get_estimated_cost(self, spot: bool = False) -> float:
        """Returns estimate cost of a single machine in the group.

        This method is an overlay of the more general method, but
        it verifies if the cost has already been estimated and returns
        it immediately if it has.
        """
        if self._estimated_cost is not None:
            return self._estimated_cost

        self._estimated_cost = inductiva.resources.estimate_machine_cost(
            self.machine_type,
            spot,
        )

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
        logging.info("Machine Type: %s", self.machine_type)
        logging.info("> Disk size: %s GB", self.disk_size_gb)
