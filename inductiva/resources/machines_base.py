"""Base class for machine groups."""
from collections import defaultdict, namedtuple
import time
import enum
import json
from typing import Optional, Union
from abc import abstractmethod
import datetime

import logging

import inductiva
import inductiva.client.models
from inductiva import api, users
from inductiva.utils import format_utils
from inductiva.client.apis.tags import compute_api

from inductiva.resources import machine_types

VCPUCount = namedtuple("VCPUCount", ["total", "per_machine"])


class ResourceType(enum.Enum):
    """Enum to represent the type of machine to be launched."""

    STANDARD = "standard"
    MPI = "mpi"


class BaseMachineGroup:
    """Base class to manage Google Cloud resources."""

    def __init__(
        self,
        machine_type: str,
        provider: Union[machine_types.ProviderType, str] = "GCP",
        data_disk_gb: int = 10,
        max_idle_time: Optional[datetime.timedelta] = None,
        auto_terminate_ts: Optional[datetime.datetime] = None,
        register: bool = True,
    ) -> None:
        """Create a BaseMachineGroup object.

        Args:
            machine_type: The type of GC machine to launch. Ex: "e2-standard-4".
              Check https://cloud.google.com/compute/docs/machine-resource for
              more information about machine types.
            provider: The cloud provider of the machine group.
            data_disk_gb: The size of the disk for user data (in GB).
            max_idle_time: Time without executing any task, after which the
              resource will be terminated.
            auto_terminate_ts: Moment in which the resource will be
              automatically terminated.
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
        self._max_idle_time = max_idle_time
        self._auto_terminate_ts = auto_terminate_ts

        # Set the API configuration that carries the information from the client
        # to the backend.
        self._api = compute_api.ComputeApi(api.get_client())
        self._estimated_cost = None

    @property
    def id(self):
        return self._id

    @property
    def n_vcpus(self):
        """Returns the number of vCPUs available in the resource.

        Returns a tuple with the total number of vCPUs and the number of vCPUs
        per machine. For a machine group with 2 machines, each with 4 vCPUs,
        this will return (8, 4).
        In a case of an Elastic machine group, the total number of vCPUs is
        the maximum number of vCPUs that can be used at the same time.
        """
        max_vcpus = int(self.quota_usage["max_vcpus"])
        max_instances = int(self.quota_usage["max_instances"])

        return VCPUCount(max_vcpus, max_vcpus // max_instances)

    @property
    def name(self):
        return self._name

    @property
    def max_idle_time(self) -> datetime.timedelta:
        return self._max_idle_time

    @property
    def auto_terminate_ts(self):
        return self._auto_terminate_ts

    @property
    def idle_time(self) -> datetime.timedelta:
        """
        Resource idle time in seconds.
        """
        return self._idle_seconds

    @staticmethod
    def _timedelta_to_seconds(
            value: Optional[datetime.timedelta] = None) -> Optional[float]:
        """Converts a timedelta object to seconds."""
        return value.total_seconds() if value is not None else None

    @staticmethod
    def _seconds_to_timedelta(
            value: Optional[float] = None) -> Optional[datetime.timedelta]:
        """Converts seconds to a timedelta object."""
        return datetime.timedelta(
            seconds=float(value)) if value is not None else None

    @staticmethod
    def _convert_auto_terminate_ts(
            timestamp: Optional[datetime.datetime]) -> Optional[str]:
        """Converts a datetime object to ISO format. Additionally, it checks
        if the datetime is in the future and if it is timezone aware."""
        if timestamp is not None:
            if (timestamp.tzinfo is None or
                    timestamp.tzinfo.utcoffset(timestamp) is None):
                raise ValueError("auto_terminate_ts must be timezone aware.")

            now_ts = datetime.datetime.now(timestamp.tzinfo)
            if timestamp < now_ts:
                raise ValueError("auto_terminate_ts must be in the future.")
            return timestamp.isoformat()

        return None

    @staticmethod
    def _iso_to_datetime(
            timestamp: Optional[str]) -> Optional[datetime.datetime]:
        """Converts an ISO format string back to a datetime object. It ensures
        the datetime is timezone aware."""
        if timestamp is not None:
            dt = datetime.datetime.fromisoformat(str(timestamp))
            if dt.tzinfo is None or dt.tzinfo.utcoffset(dt) is None:
                raise ValueError("The datetime string must be timezone aware.")
            return dt
        return None

    def _register_machine_group(self, **kwargs):
        """Register machine group configuration in API.

        Returns:
            The unique ID and name identifying the machine on the API."""

        instance_group_config = inductiva.client.models.VMGroupConfig(
            machine_type=self.machine_type,
            provider_id=self.provider,
            disk_size_gb=self.data_disk_gb,
            max_idle_time=self._timedelta_to_seconds(self.max_idle_time),
            auto_terminate_ts=self._convert_auto_terminate_ts(
                self.auto_terminate_ts),
            **kwargs,
        )

        resp = self._api.register_vm_group(
            body=instance_group_config,
            skip_deserialization=True,
        )
        body = json.loads(resp.response.data)

        self._id = body["id"]
        self._name = body["name"]
        self.quota_usage = body.get("quota_usage") or {}
        self.register = False
        # Lifecycle configuration parameters are updated with default values
        # from the API response if they were not provided by the user
        self._max_idle_time = self._seconds_to_timedelta(
            body.get("max_idle_time"))
        self._idle_seconds = self._seconds_to_timedelta(
            body.get("idle_seconds"))
        self._auto_terminate_ts = self._iso_to_datetime(
            body.get("auto_terminate_ts"))
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
        machine_group.quota_usage = resp.get("quota_usage") or {}
        machine_group._max_idle_time = cls._seconds_to_timedelta(
            resp.get("max_idle_time"))
        machine_group._idle_seconds = cls._seconds_to_timedelta(
            resp.get("idle_seconds"))
        machine_group._auto_terminate_ts = cls._iso_to_datetime(
            resp.get("auto_terminate_ts"))

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
        self._api.start_vm_group(body=request_body)
        creation_time = format_utils.seconds_formatter(time.time() - start_time)
        self._started = True
        logging.info("%s successfully started in %s.", self, creation_time)
        logging.info("")
        logging.info("The machine group is using the following quotas:")
        self.log_quota_usage("used by resource")
        return True

    def terminate(self, verbose: bool = True, **kwargs):
        """Terminates a machine group."""
        if not self._started or self.id is None or self.name is None:
            logging.warning(
                "Attempting to terminate an unstarted machine group.")
            return

        try:
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
            if verbose:
                logging.info("Successfully requested termination of %s.",
                             repr(self))
                logging.info("Termination of the machine group "
                             "freed the following quotas:")
                self.log_quota_usage("freed by resource")
            return True

        except inductiva.client.ApiException as api_exception:
            raise api_exception

    def _get_estimated_cost(self, spot: bool = True) -> float:
        """Returns estimate cost of a single machine in the group.

        This method is an overlay of the more general method, but
        it verifies if the cost has already been estimated and returns
        it immediately if it has.
        """
        if self.provider in (
                machine_types.ProviderType.ICE,
                machine_types.ProviderType.LOCAL,
        ):
            return 0

        self._estimated_cost = inductiva.resources.estimate_machine_cost(
            self.machine_type,
            spot,
        )

        return self._estimated_cost

    def log_quota_usage(self, resource_usage_header: str):
        quotas = users.get_quotas()
        table = defaultdict(list)
        emph_formatter = format_utils.get_ansi_formatter()

        header_formatters = [
            lambda x: emph_formatter(x.upper(), format_utils.Emphasis.BOLD)
        ]

        for name, value in self.quota_usage.items():
            max_allowed = quotas.get(name, {}).get("max_allowed", "n/a")
            full_name = quotas.get(name, {}).get("label", "n/a")
            in_use = quotas.get(name, {}).get("in_use", "n/a")

            table[""].append(full_name)
            table[resource_usage_header].append(value)
            table["current usage"].append(in_use)
            table["max allowed"].append(max_allowed)

        table_str = format_utils.get_tabular_str(
            table,
            header_formatters=header_formatters,
        )

        logging.info(table_str)

    def _log_estimated_spot_vm_savings(self) -> None:
        if self.provider in (
                machine_types.ProviderType.ICE,
                machine_types.ProviderType.LOCAL,
        ):
            return

        spot_cost = self._get_estimated_cost(True)
        non_spot_cost = self._get_estimated_cost(False)
        spot_times_cheaper = round(non_spot_cost / spot_cost, 2)

        is_spot = getattr(self, "spot", True)

        if not is_spot:
            logging.info(
                "\t· The same machine group with spot machines would cost "
                "%.1fx less. Specify "
                "`spot=True` in the constructor to use spot machines.",
                spot_times_cheaper,
            )
        else:
            logging.info(
                "\t· You are spending %.1fx less "
                "by using spot machines.",
                spot_times_cheaper,
            )
        logging.info("")

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

        logging.info("\t· Name:                       %s", self.name)
        logging.info("\t· Machine Type:               %s", self.machine_type)
        logging.info("\t· Data disk size:             %s GB", self.data_disk_gb)

        # Log max idle time
        value_str = format_utils.timedelta_formatter(
            self.max_idle_time) if self.max_idle_time is not None else "N/A"
        logging.info("\t· Maximum idle time:          %s", value_str)

        # Log auto terminate timestamp
        value_str = self.auto_terminate_ts.strftime(
            "%Y/%m/%d %H:%M:%S"
        ) if self.auto_terminate_ts is not None else "N/A"
        logging.info("\t· Auto terminate timestamp:   %s", value_str)
