"""Base class for machine groups."""
from collections import defaultdict, namedtuple
from typing import Optional, Union
from abc import abstractmethod
from abc import ABC
import datetime
import time
import enum
import json
import math

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


class BaseMachineGroup(ABC):
    """Base class to manage Google Cloud resources."""

    QUOTAS_EXCEEDED_SLEEP_SECONDS = 60

    def __init__(
        self,
        machine_type: str,
        provider: Union[machine_types.ProviderType, str] = "GCP",
        threads_per_core: int = 2,
        data_disk_gb: int = 10,
        auto_resize_disk_max_gb: int = 500,
        max_idle_time: Optional[Union[datetime.timedelta, int]] = None,
        auto_terminate_ts: Optional[datetime.datetime] = None,
        register: bool = True,
    ) -> None:
        """Create a BaseMachineGroup object.

        Args:
            machine_type: The type of GC machine to launch. Ex: "e2-standard-4".
              Check https://cloud.google.com/compute/docs/machine-resource for
              more information about machine types.
            provider: The cloud provider of the machine group.
            threads_per_core: The number of threads per core (1 or 2).
            data_disk_gb: The size of the disk for user data (in GB).
            auto_resize_disk_max_gb: The maximum size in GB that the hard disk
                of the cloud VM can reach. If set, the disk will be
                automatically resized, during the execution of a task, when the
                free space falls below a certain threshold. This mechanism helps
                prevent "out of space" errors, that can occur when a task
                generates a quantity of output files that exceeds the size of
                the local storage. Increasing disk size during task execution
                increases the cost of local storage associated with the VM,
                therefore the user must set an upper limit to the disk size, to
                prevent uncontrolled costs. Once that limit is reached, the disk
                is no longer automatically resized, and if the task continues to
                output files, it will fail.
            max_idle_time: Time without executing any task, after which the
              resource will be terminated. Can be an exact timedelta or an int
                representing the number of minutes.
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
        self._free_space_threshold_gb = 5
        self._size_increment_gb = 10

        if data_disk_gb <= 0:
            raise ValueError("`data_disk_gb` must be positive.")

        if auto_resize_disk_max_gb is not None:
            if not isinstance(auto_resize_disk_max_gb,
                              int) or auto_resize_disk_max_gb <= 0:
                raise ValueError(
                    "`auto_resize_disk_max_gb` must be a positive integer.")

            if auto_resize_disk_max_gb < data_disk_gb + self._size_increment_gb:
                raise ValueError("`auto_resize_disk_max_gb` must be greater "
                                 "than or equal to `data_disk_gb + "
                                 f"{self._size_increment_gb}GB`.")

        if threads_per_core not in [1, 2]:
            raise ValueError("`threads_per_core` must be either 1 or 2.")

        self.machine_type = machine_type
        self.provider = provider.value
        self.threads_per_core = threads_per_core
        self.data_disk_gb = data_disk_gb
        self.auto_resize_disk_max_gb = auto_resize_disk_max_gb
        self._id = None
        self._name = None
        self.create_time = None
        self._started = False
        self.register = register
        #Number of active machines at the time of
        #the request machine_groups.get()
        self._active_machines = 0
        self.num_machines = 0
        self._auto_terminate_ts = auto_terminate_ts
        self._custom_vm_image = None

        # Set the API configuration that carries the information from the client
        # to the backend.
        self._api = compute_api.ComputeApi(api.get_client())
        self._estimated_cost = None
        self._max_idle_time = max_idle_time

        if isinstance(max_idle_time, int):
            if max_idle_time <= 0:
                raise ValueError("`max_idle_time` must be positive.")
            self._max_idle_time = datetime.timedelta(minutes=max_idle_time)

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

        # if threads per core is 2
        cores_per_machine = max_vcpus // max_instances

        if self.threads_per_core == 1:
            cores_per_machine //= 2

        return VCPUCount(cores_per_machine * max_instances, cores_per_machine)

    @property
    def available_vcpus(self):
        """Returns the maximum number of vCPUs that can be used on a task.

        On a machine group with 2 machines, each with 4 vCPUs, this will return
        4.
        On an elastic machine group, this will also return 4.
        On an MPI cluster this will return the total number of vcpus because
        we can run on the total number of vcpus.
        """
        return self.n_vcpus.per_machine

    @property
    def name(self):
        return self._name

    @property
    def started(self):
        return self._started

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
        return datetime.timedelta(seconds=float(value)) if value else None

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

    def _dynamic_disk_resize_config(self):
        if self.auto_resize_disk_max_gb is None:
            return None

        return {
            "free_space_threshold_gb": self._free_space_threshold_gb,
            "size_increment_gb": self._size_increment_gb,
            "max_disk_size_gb": self.auto_resize_disk_max_gb
        }

    @staticmethod
    def _iso_to_datetime(
            timestamp: Optional[str]) -> Optional[datetime.datetime]:
        """Converts an ISO format string back to a datetime object. It ensures
        the datetime is timezone aware."""
        if timestamp:
            dt = datetime.datetime.fromisoformat(str(timestamp))
            if dt.year == 9999:
                return None

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
            threads_per_core=self.threads_per_core,
            disk_size_gb=self.data_disk_gb,
            max_idle_time=self._timedelta_to_seconds(self.max_idle_time),
            auto_terminate_ts=self._convert_auto_terminate_ts(
                self.auto_terminate_ts),
            dynamic_disk_resize_config=self._dynamic_disk_resize_config(),
            custom_vm_image=self._custom_vm_image,
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
        self.total_ram_gb = body.get("total_ram_gb")
        self._cost_per_hour = body.get("cost_per_hour")

        dynamic_disk_resize_config = body.get(
            "dynamic_disk_resize_config") or {}
        self.auto_resize_disk_max_gb = dynamic_disk_resize_config.get(
            "max_disk_size_gb")
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

    def can_start_resource(self) -> bool:
        """Check if the resource can be started.

        This method checks if the resource can be started by checking
        the available quotas and resource usage.

        returns:
            bool: True if the resource can be started, False otherwise.
        """
        quotas = users.get_quotas()

        cost_in_use = quotas["max_price_hour"]["in_use"]
        cost_max = math.inf if quotas["max_price_hour"][
            "max_allowed"] is None else quotas["max_price_hour"]["max_allowed"]
        estimated_cost = cost_in_use + self.estimate_cloud_cost(verbose=False)
        is_cost_ok = estimated_cost <= cost_max

        vcpu_in_use = quotas["max_vcpus"]["in_use"]
        vcpu_max = math.inf if quotas["max_vcpus"][
            "max_allowed"] is None else quotas["max_vcpus"]["max_allowed"]
        current_vcpu = self.n_vcpus.total
        estimated_vcpu_usage = vcpu_in_use + current_vcpu
        is_vcpu_ok = estimated_vcpu_usage <= vcpu_max

        machines_in_use = quotas["max_instances"]["in_use"]
        machines_max = math.inf if quotas["max_instances"][
            "max_allowed"] is None else quotas["max_instances"]["max_allowed"]

        estimated_machine_usage = machines_in_use + self.num_machines
        is_instance_ok = estimated_machine_usage <= machines_max

        return is_cost_ok and is_vcpu_ok and is_instance_ok

    def start(self, wait_for_quotas: bool = False, **kwargs):
        """Starts a machine group.

        Args:
            wait_for_quotas: If True, the method will wait for quotas to
              become available before starting the resource.
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
                threads_per_core=self.threads_per_core,
                disk_size_gb=self.data_disk_gb,
                **kwargs,
            )

        logging.info(
            "Starting %s. This may take a few minutes.\n"
            "Note that stopping this local process will not interrupt "
            "the creation of the machine group. Please wait...", repr(self))
        start_time = time.time()

        if wait_for_quotas:
            if not self.can_start_resource():
                print("This machine will exceed the current quotas.\n"
                      "Will wait for quotas to become available.")
            while not self.can_start_resource():
                time.sleep(self.QUOTAS_EXCEEDED_SLEEP_SECONDS)

        self._api.start_vm_group(body=request_body)
        creation_time = format_utils.seconds_formatter(time.time() - start_time)
        self._started = True
        quota_usage_table_str = self.quota_usage_table_str("used by resource")
        logging.info(
            "%s successfully started in %s.\n\n"
            "The machine group is using the following quotas:\n"
            "%s", self, creation_time, quota_usage_table_str)
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
                    threads_per_core=self.threads_per_core,
                    disk_size_gb=self.data_disk_gb,
                    **kwargs,
                )

            self._api.delete_vm_group(body=request_body)
            if verbose:
                logging.info("Successfully requested termination of %s.",
                             repr(self))
                logging.info("Termination of the machine group "
                             "freed the following quotas:")
                logging.info(self.quota_usage_table_str("freed by resource"))
            return True

        except inductiva.client.ApiException as api_exception:
            raise api_exception

    def _get_estimated_cost(self, spot: bool = True) -> float:
        """Returns estimate cost of a single machine in the group.

        This method is an overlay of the more general method, but
        it verifies if the cost has already been estimated and returns
        it immediately if it has.
        """
        if self.provider in (machine_types.ProviderType.LOCAL,):
            return 0

        self._estimated_cost = inductiva.resources.estimate_machine_cost(
            self.machine_type,
            spot,
        )

        return self._estimated_cost

    def quota_usage_table_str(self, resource_usage_header: str) -> str:
        quotas = users.get_quotas()
        table = defaultdict(list)
        emph_formatter = format_utils.get_ansi_formatter()

        header_formatters = [
            lambda x: emph_formatter(x.upper(), format_utils.Emphasis.BOLD)
        ]

        def _format_float(x):
            return f"{x:.3f}" if isinstance(x, float) else str(x)

        formatters = {
            column: [_format_float] for column in
            [resource_usage_header, "current usage", "max allowed"]
        }

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
            formatters=formatters,
            header_formatters=header_formatters,
        )

        return table_str

    def _log_estimated_spot_vm_savings(self) -> None:
        if self.provider in (machine_types.ProviderType.LOCAL,):
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

        if self.auto_resize_disk_max_gb:
            logging.info("\t· Auto resize disk max size:  %s GB",
                         self.auto_resize_disk_max_gb)

        logging.info("\t· Total memory (RAM):         %s GB", self.total_ram_gb)

        # Log max idle time
        value_str = format_utils.timedelta_formatter(
            self.max_idle_time) if self.max_idle_time is not None else "N/A"
        logging.info("\t· Maximum idle time:          %s", value_str)

        # Log auto terminate timestamp
        value_str = self.auto_terminate_ts.strftime(
            "%Y/%m/%d %H:%M:%S"
        ) if self.auto_terminate_ts is not None else "N/A"
        logging.info("\t· Auto terminate timestamp:   %s", value_str)

    def estimate_cloud_cost(self, verbose: bool = True):
        """Estimates a cost per hour of min and max machines in US dollars.

        these are the estimted costs of having minimum and the
        maximum number of machines up in the cloud. The final cost will vary
        depending on the total usage of the machines."""

        min_cost_per_hour = self._cost_per_hour.get("min")
        max_cost_per_hour = self._cost_per_hour.get("max")

        if not verbose:
            return max_cost_per_hour

        if min_cost_per_hour == max_cost_per_hour:
            logging.info(
                "\t· Estimated cloud cost of machine group: %.3f $/h",
                max_cost_per_hour,
            )
        else:
            min_reason = self._cost_per_hour.get("min_reason").rstrip(".")
            max_reason = self._cost_per_hour.get("max_reason").rstrip(".")

            logging.info("\t· Estimated cloud cost of machine group:")
            logging.info("\t\t· Minimum: %.3f $/h (%s)", min_cost_per_hour,
                         min_reason)
            logging.info("\t\t· Maximum: %.3f $/h (%s)", max_cost_per_hour,
                         max_reason)

        self._log_estimated_spot_vm_savings()

        return max_cost_per_hour
