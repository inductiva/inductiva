"""Base class for machine groups."""
import time
import enum
from abc import abstractmethod

import logging

import inductiva
import inductiva.client.models
from inductiva import api
from inductiva.utils import format_utils
from inductiva.client.apis.tags import compute_api
from inductiva.client import exceptions
from inductiva import logs


class ResourceType(enum.Enum):
    """Enum to represent the type of machine to be launched."""

    STANDARD = "standard"
    MPI = "mpi"


class BaseMachineGroup:
    """Base class to manage Google Cloud resources."""

    def __init__(self,
                 machine_type: str,
                 data_disk_gb: int = 10,
                 register: bool = True) -> None:
        """Create a BaseMachineGroup object.

        Args:
            machine_type: The type of GC machine to launch. Ex: "e2-standard-4".
              Check https://cloud.google.com/compute/docs/machine-resource for
              more information about machine types.
            spot: Whether to use spot machines.
            data_disk_gb: The size of the disk for user data (in GB).
            register: Bool that indicates if a machine group should be register
                or if it was already registered. If set to False by users on
                initialization, then, the machine group will not be able to be
                started. This serves has an helper argument for retrieving
                already registered machine groups that can be started, for
                example, when retrieving with the `machines_groups.get` method.
                Users should not set this argument in anyway.
        """
        if machine_type not in inductiva.resources.list_available_machines():
            raise ValueError("Machine type not supported")

        if data_disk_gb < 0 or data_disk_gb > 100:
            raise ValueError(
                "`data_disk_gb` must be a positive value smaller than 100 GB")

        self.machine_type = machine_type
        self.data_disk_gb = data_disk_gb
        self._true_disk_size_gb =\
            data_disk_gb + inductiva.constants.BASE_MACHINE_DISK_SIZE_GB
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
            disk_size_gb=self._true_disk_size_gb,
            **kwargs,
        )

        try:
            resp = self._api.register_vm_group(body=instance_group_config)
        except (exceptions.ApiValueError, exceptions.ApiException) as e:
            logs.log_and_exit(
                logging.getLogger(),
                logging.ERROR,
                "Registering machine group failed with exception %s",
                e,
                exc_info=e)

        self._id = resp.body["id"]
        self._name = resp.body["name"]
        self.register = False
        self._log_machine_group_info()

    @abstractmethod
    def __repr__(self):
        pass

    @classmethod
    def from_api_response(cls, resp: dict):
        """Creates a MachineGroup object from an API response."""

        machine_group = cls(
            machine_type=resp["machine_type"],
            data_disk_gb=resp["disk_size_gb"] -
            inductiva.constants.BASE_MACHINE_DISK_SIZE_GB,
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
                disk_size_gb=self._true_disk_size_gb,
                **kwargs,
            )
        logging.info("Starting %s. "
                     "This may take a few minutes.", repr(self))
        logging.info("Note that stopping this local process will not interrupt "
                     "the creation of the machine group. Please wait...")
        start_time = time.time()
        try:
            self._api.start_vm_group(body=request_body)
        except inductiva.client.ApiException as e:
            logs.log_and_exit(logging.getLogger(),
                              logging.ERROR,
                              "Starting machine group failed with exception %s",
                              e,
                              exc_info=e)
        creation_time = format_utils.seconds_formatter(time.time() - start_time)
        self._started = True
        logging.info("%s successfully started in %s.", self, creation_time)

    def terminate(self, **kwargs):
        """Terminates a machine group."""
        if not self._started or self.id is None or self.name is None:
            logging.info("Attempting to terminate an unstarted machine group.")
            return

        try:
            logging.info("Terminating %s. This may take a few minutes.",
                         repr(self))
            start_time = time.time()

            request_body = \
                inductiva.client.models.GCPVMGroup(
                    id=self.id,
                    name=self.name,
                    machine_type=self.machine_type,
                    disk_size_gb=self._true_disk_size_gb,
                    **kwargs,
                )

            self._api.delete_vm_group(body=request_body)
            termination_time = format_utils.seconds_formatter(time.time() -
                                                              start_time)
            logging.info("%s successfully terminated in %s.", self,
                         termination_time)

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
                "Attempting to get the status of an unregistered machine group."
            )
            return

        response = self._api.get_group_status({"name": self.name})
        if response.body == "notFound":
            logging.info("Machine group does not exist: %s.", self.name)
        return response.body

    def _log_machine_group_info(self):
        """Logs the machine group info."""

        logging.info("> Name:         %s", self.name)
        logging.info("> Machine Type: %s", self.machine_type)
        logging.info("> Data disk size:    %s GB", self.data_disk_gb)
