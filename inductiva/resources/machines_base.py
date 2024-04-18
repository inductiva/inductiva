"""Base class for machine groups."""
import time
import enum
from typing import Union
from abc import abstractmethod
import datetime

import logging

import inductiva
import inductiva.client.models
from inductiva import api
from inductiva.utils import format_utils
from inductiva.client.apis.tags import compute_api
from inductiva.client import exceptions
from inductiva import logs

from inductiva.resources import machine_types


class ResourceType(enum.Enum):
    """Enum to represent the type of machine to be launched."""

    STANDARD = "standard"
    MPI = "mpi"


class BaseMachineGroup:
    """Base class to manage Google Cloud resources."""

    def __init__(self,
                 machine_type: str,
                 provider: Union[machine_types.ProviderType, str] = "GCP",
                 data_disk_gb: int = 10,
                 register: bool = True) -> None:
        """Create a BaseMachineGroup object.

        Args:
            machine_type: The type of GC machine to launch. Ex: "e2-standard-4".
              Check https://cloud.google.com/compute/docs/machine-resource for
              more information about machine types.
            data_disk_gb: The size of the disk for user data (in GB).
            register: Bool that indicates if a machine group should be register
                or if it was already registered. If set to False by users on
                initialization, then, the machine group will not be able to be
                started. This serves has an helper argument for retrieving
                already registered machine groups that can be started, for
                example, when retrieving with the `machines_groups.get` method.
                Users should not set this argument in anyway.
        """

        provider = machine_types.ProviderType(provider)
        self.provider = provider.value

        if machine_type not in machine_types.list_available_machines(
                self.provider):
            raise ValueError(f"Machine type not supported in {self.provider}")

        if data_disk_gb <= 0:
            raise ValueError("`data_disk_gb` must be positive.")

        self.machine_type = machine_type
        self.provider = provider.value
        self.data_disk_gb = data_disk_gb
        self._id = None
        self._name = None
        self.create_time = None
        self._started = False
        self.register = register
        #Number of active machines at the time of
        #the request machine_groups.get()
        self._active_machines = 0
        self.num_machines = 0

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

        instance_group_config = inductiva.client.models.VMGroupConfig(
            machine_type=self.machine_type,
            provider_id=self.provider,
            disk_size_gb=self.data_disk_gb,
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

    def active_machines_to_str(self) -> str:
        """Return the number of machines currently running.
        """
        return f"{self._active_machines}/{self.num_machines}"

    @classmethod
    def from_api_response(cls, resp: dict):
        """Creates a MachineGroup object from an API response."""

        machine_group = cls(
            machine_type=resp["machine_type"],
            data_disk_gb=resp["disk_size_gb"],
            register=False,
        )
        machine_group._id = resp["id"]
        machine_group.provider = resp["provider_id"]
        machine_group._name = resp["name"]
        machine_group.create_time = resp["creation_timestamp"]
        machine_group._started = bool(resp["started"])

        return machine_group

    def update_termination_timers(self,
                                  max_idle_time: datetime.timedelta = None,
                                  auto_terminate_ts: datetime.datetime = None):
        """Update the termination timers of a machine group.

        Args:
            max_idle_time (timedelta): Max idle time, i.e. time without
                executing any task, after which the resource will be terminated.
            auto_terminate_ts (datetime): Moment in which the resource will
                be automatically terminated, irrespectively of the existence of
                tasks yet to be executed by the resource.
        """

        # Convert max_idle_time from minutes to seconds
        max_idle_time = max_idle_time.total_seconds() if max_idle_time else None

        # Convert auto_terminate_ts to ISO format
        if auto_terminate_ts is not None:
            now_ts = datetime.datetime.now(auto_terminate_ts.tzinfo)
            if auto_terminate_ts < now_ts:
                raise ValueError("auto_terminate_ts must be in the future.")
            auto_terminate_ts = auto_terminate_ts.isoformat()

        update_body = inductiva.client.models.VMGroupLifecycleConfig(
            max_idle_time=max_idle_time, auto_terminate_ts=auto_terminate_ts)

        if max_idle_time is not None or auto_terminate_ts is not None:
            try:
                self._api.update_vm_group_config(
                    path_params={"mg_id": self._id}, body=update_body)
            except inductiva.client.ApiException as e:
                logs.log_and_exit(
                    logging.getLogger(),
                    logging.ERROR,
                    "Setting termination timers for machine group failed " \
                    "with exception %s.",
                    e,
                    exc_info=e)

    def start(self,
              max_idle_time: datetime.timedelta = None,
              auto_terminate_ts: datetime.datetime = None,
              **kwargs):
        """Starts a machine group.

        Args:
            max_idle_time (timedelta): Max idle time, i.e. time without
                executing any task, after which the resource will be terminated.
            auto_terminate_ts (datetime): Moment in which the resource will
                be automatically terminated, irrespectively of the existence of
                tasks yet to be executed by the resource.
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
            inductiva.client.models.VMGroupConfig(
                id=self.id,
                name=self.name,
                machine_type=self.machine_type,
                provider_id=self.provider,
                disk_size_gb=self.data_disk_gb,
                **kwargs,
            )

        logging.info("Starting %s. "
                     "This may take a few minutes.", repr(self))
        logging.info("Note that stopping this local process will not interrupt "
                     "the creation of the machine group. Please wait...")
        start_time = time.time()

        self.update_termination_timers(max_idle_time, auto_terminate_ts)

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
                inductiva.client.models.VMGroupConfig(
                    id=self.id,
                    name=self.name,
                    machine_type=self.machine_type,
                    provider_id=self.provider,
                    disk_size_gb=self.data_disk_gb,
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
        if self.provider == "ICE":
            return 0

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
