"""Base class for machine groups."""
from collections import defaultdict, namedtuple
from dataclasses import dataclass
from typing import Optional, Union
from abc import ABC, abstractmethod
import datetime
import time
import enum
import json
import math

import logging

import inductiva
import inductiva.client.models
from inductiva import api, users
from inductiva.resources.utils import ProviderType
from inductiva.utils import format_utils
from inductiva.client.apis.tags import compute_api

VCPUCount = namedtuple("VCPUCount", ["total", "per_machine"])


class ResourceType(enum.Enum):
    """Enum to represent the type of machine to be launched."""

    STANDARD = "standard"
    MPI = "mpi"


@dataclass(repr=False)
class BaseMachineGroup(ABC):
    """Base class to manage Google Cloud resources.

    Args:
        machine_type: The type of GC machine to launch. Ex: "e2-standard-4".
          Check https://cloud.google.com/compute/docs/machine-resource for
          more information about machine types.
        zone: The zone where the machines will be launched.
        provider: The cloud provider of the machine group.
        threads_per_core: The number of threads per core (1 or 2).
        data_disk_gb: The size of the disk for user data (in GB).
        max_idle_time: Time without executing any task, after which the
          resource will be terminated. Can be an exact timedelta or an int
            representing the number of minutes.
        auto_terminate_ts: Moment in which the resource will be
          automatically terminated.
        auto_terminate_minutes: Duration, in minutes, the MPICluster will be
                kept alive. After auto_terminate_minutes minutes the machine
                will be terminated. This time will start counting after calling
                this method.

    :meta private:
    """
    # Constructor arguments
    machine_type: str
    zone: Optional[str] = "europe-west1-b"
    provider: Union[ProviderType, str] = "GCP"
    threads_per_core: int = 2
    data_disk_gb: int = 10
    max_idle_time: Union[datetime.timedelta, int] = 3
    auto_terminate_ts: Optional[datetime.datetime] = None
    auto_terminate_minutes: Optional[int] = None

    create_time = None
    num_machines = 0
    quota_usage = {}
    allow_auto_start = True

    # Internal attributes
    _free_space_threshold_gb = 5
    _id = None
    _name = None
    _started = False
    #Number of active machines at the time of
    #the request machine_groups.get()
    _active_machines = 0
    _custom_vm_image = None
    _estimated_cost = None
    _idle_seconds = None
    _cost_per_hour = {}
    _total_ram_gb = None
    _cpu_info = {}
    _gpu_info = {}

    QUOTAS_EXCEEDED_SLEEP_SECONDS = 60

    def __post_init__(self):
        """Validate inputs and initialize additional attributes after
        dataclass initialization."""
        provider = ProviderType(self.provider)
        self.provider = provider.value

        # Set the API configuration that carries the information from the client
        # to the backend.
        self._api = compute_api.ComputeApi(api.get_client())

        self._validate_inputs()

    def _validate_inputs(self):
        """Validate initialization inputs."""
        if not isinstance(self.data_disk_gb, int):
            raise ValueError("`data_disk_gb` must be an integer.")

        if not 10 <= self.data_disk_gb <= 65536:
            raise ValueError("`data_disk_gb` must be between 10 and 65536.")

        if self.auto_resize_disk_max_gb is not None:
            if not isinstance(self.auto_resize_disk_max_gb,
                              int) or self.auto_resize_disk_max_gb <= 0:
                raise ValueError(
                    "`auto_resize_disk_max_gb` must be a positive integer.")

            if self.auto_resize_disk_max_gb < self.data_disk_gb:
                raise ValueError(
                    "`auto_resize_disk_max_gb` must be greater than \
                    or equal to `data_disk_gb GB`.")

        if self.threads_per_core not in [1, 2]:
            raise ValueError("`threads_per_core` must be either 1 or 2.")

        if isinstance(self.max_idle_time, int):
            if self.max_idle_time <= 0:
                raise ValueError("`max_idle_time` must be positive.")
            self._max_idle_time = datetime.timedelta(minutes=self.max_idle_time)

        if self.auto_terminate_ts is not None:
            logging.warning("You are using `auto_terminate_ts`. This argument"
                            "will be deprecated in the future. Please use"
                            "`auto_terminate_minutes` instead.")

        if isinstance(self.auto_terminate_minutes, int):
            time_delta_minutes = datetime.timedelta(
                minutes=self.auto_terminate_minutes)
            self._auto_terminate_ts = datetime.datetime.now(
                tz=datetime.timezone.utc) + time_delta_minutes

    def has_gpu(self) -> bool:
        """Check if the machine group has a GPU."""
        if self._gpu_info is None:
            return False
        return self._gpu_info.get("gpu_count") > 0

    @property
    def id(self):
        return self._id

    @property
    @abstractmethod
    def n_vcpus(self):
        """Returns the number of vCPUs available in the resource.

        Returns a tuple with the total number of vCPUs and the number of vCPUs
        per machine. For a machine group with 2 machines, each with 4 vCPUs,
        this will return (8, 4).
        In a case of an Elastic machine group, the total number of vCPUs is
        the maximum number of vCPUs that can be used at the same time.
        """
        pass

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
    def idle_time(self) -> datetime.timedelta:
        """
        Resource idle time in seconds.
        """
        return self._idle_seconds

    @property
    def total_ram_gb(self):
        return self._total_ram_gb

    @abstractmethod
    def short_name(self) -> str:
        pass

    @staticmethod
    def _timedelta_to_seconds(value: Union[datetime.timedelta, int]) -> float:
        """Converts a timedelta object to seconds."""
        if isinstance(value, int):
            return value * 60
        return value.total_seconds()

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

    def _update_attributes_from_response(self, resp: dict):
        """Update machine group attributes with values from the API response."""
        self._id = resp["id"]
        self._name = resp["name"]
        self.quota_usage = resp.get("quota_usage") or {}
        # Lifecycle configuration parameters are updated with default values
        # from the API response if they were not provided by the user
        self.max_idle_time = self._seconds_to_timedelta(
            resp.get("max_idle_time"))
        self._idle_seconds = self._seconds_to_timedelta(
            resp.get("idle_seconds"))
        self.auto_terminate_ts = self._iso_to_datetime(
            resp.get("auto_terminate_ts"))
        self._total_ram_gb = resp.get("total_ram_gb")
        self._cost_per_hour = resp.get("cost_per_hour")
        self._cpu_info = resp.get("cpu_info")
        self._gpu_info = resp.get("gpu_info")
        self.zone = resp.get("zone")
        dynamic_disk_resize_config = resp.get(
            "dynamic_disk_resize_config") or {}
        self.auto_resize_disk_max_gb = dynamic_disk_resize_config.get(
            "max_disk_size_gb")

    def _register_machine_group(self, **kwargs):
        """Register machine group configuration in API.

        Returns:
            The unique ID and name identifying the machine on the API."""
        logging.info("■ Registering %s configurations:", self.short_name())

        instance_group_config = inductiva.client.models.RegisterVMGroupRequest(
            machine_type=self.machine_type,
            provider_id=self.provider,
            threads_per_core=self.threads_per_core,
            disk_size_gb=self.data_disk_gb,
            max_idle_time=self._timedelta_to_seconds(self.max_idle_time),
            auto_terminate_ts=self._convert_auto_terminate_ts(
                self.auto_terminate_ts),
            dynamic_disk_resize_config=self._dynamic_disk_resize_config(),
            custom_vm_image=self._custom_vm_image,
            zone=self.zone,
            **kwargs,
        )

        resp = self._api.register_vm_group(
            body=instance_group_config,
            skip_deserialization=True,
        )
        body = json.loads(resp.response.data)

        self._update_attributes_from_response(body)

        self._log_machine_group_info()

    def __repr__(self):
        class_name = self.__class__.__name__
        return f"{class_name}(name=\"{self.name}\")"

    def active_machines_to_str(self) -> str:
        """Return the number of machines currently running.
        """
        return f"{self._active_machines}/{self.num_machines}"

    @classmethod
    def from_api_response(cls, resp: dict):
        """Creates a MachineGroup object from an API response."""

        # Do not call __init__ to prevent registration of the machine group
        machine_group = cls.__new__(cls)
        machine_group._api = compute_api.ComputeApi(api.get_client())
        machine_group.machine_type = resp["machine_type"]
        machine_group.data_disk_gb = resp["disk_size_gb"]
        machine_group.provider = resp["provider_id"]
        machine_group.create_time = resp["creation_timestamp"]
        machine_group._started = bool(resp["started"])
        machine_group.__dict__["machines"] = resp["machines"]
        machine_group.__dict__["_active_machines"] = int(resp["num_vms"])

        machine_group._update_attributes_from_response(resp)

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

    def start(self, wait_for_quotas: bool = False):
        """Starts a machine group.

        Args:
            wait_for_quotas: If True, the method will wait for quotas to
              become available before starting the resource."""
        if self._started:
            logging.info("Attempting to start a machine group already started.")
            return

        if self.id is None or self.name is None:
            logging.info("Attempting to start an unregistered machine group. "
                         "Make sure you have called the constructor.")
            return

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

        self._api.start_vm_group(query_params={"machine_group_id": self.id})
        creation_time = format_utils.seconds_formatter(time.time() - start_time)
        self._started = True
        quota_usage_table_str = self.quota_usage_table_str("used by resource")
        logging.info(
            "%s successfully started in %s.\n\n"
            "The machine group is using the following quotas:\n"
            "%s", self, creation_time, quota_usage_table_str)
        return True

    def terminate(self, verbose: bool = True):
        """Terminates a machine group."""
        if not self._started or self.id is None or self.name is None:
            logging.warning(
                "Attempting to terminate an unstarted machine group.")
            return

        try:
            self._api.delete_vm_group(
                query_params={"machine_group_id": self.id})
            if verbose:
                logging.info("Successfully requested termination of %s.",
                             repr(self))
                if self.provider == ProviderType.GCP:
                    logging.info("Termination of the machine group "
                                 "freed the following quotas:")
                    logging.info(
                        self.quota_usage_table_str("freed by resource"))
            return True

        except inductiva.client.ApiException as api_exception:
            raise api_exception

    def _get_estimated_cost(self, spot: bool = True) -> float:
        """Returns estimate cost of a single machine in the group.

        This method is an overlay of the more general method, but
        it verifies if the cost has already been estimated and returns
        it immediately if it has.
        """
        if self.provider in (ProviderType.LOCAL,):
            return 0

        self._estimated_cost = inductiva.resources.estimate_machine_cost(
            self.machine_type,
            spot,
            self.zone,
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
        if self.provider in (ProviderType.LOCAL,):
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

    def _log_machine_group_info(self):
        """Logs the machine group info."""

        logging.info("\t· Name:                       %s", self.name)
        logging.info("\t· Provider:                   %s", self.provider)
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


@dataclass(repr=False)
class MachineGroup(BaseMachineGroup):
    """Create a MachineGroup object.
    
    A machine group is a collection of homogenous machines with given the
    configurations that are launched in Google Cloud.
    Note: The machine group will be available only after calling 'start' method.
    The billing will start only after the machines are started.

    Args:
        machine_type: The type of GC machine to launch. Ex: "e2-standard-4".
          Check https://cloud.google.com/compute/docs/machine-resource for
          information about machine types.
        zone: The zone where the machines will be launched.
        provider: The cloud provider of the machine group.
        threads_per_core: The number of threads per core (1 or 2).
        data_disk_gb: The size of the disk for user data (in GB).
        max_idle_time: Time without executing any task, after which the
          resource will be terminated.
        auto_terminate_ts: Moment in which the resource will be
          automatically terminated.
        auto_terminate_minutes: Duration, in minutes, the MPICluster will be
                kept alive. After auto_terminate_minutes minutes the machine
                will be terminated. This time will start counting after calling
                this method.
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
        num_machines: The number of virtual machines to launch.
        spot: Whether to use spot machines.
    """
    # Constructor arguments
    auto_resize_disk_max_gb: Optional[int] = None
    num_machines: int = 1
    spot: bool = True

    # Internal attributes
    _is_elastic = False

    def __post_init__(self):
        """Validate inputs and initialize additional attributes after
        dataclass initialization."""
        super().__post_init__()

        self._register_machine_group(num_vms=self.num_machines,
                                     spot=self.spot,
                                     is_elastic=self._is_elastic)

    def _validate_inputs(self):
        super()._validate_inputs()
        if self.num_machines < 1:
            raise ValueError(
                "`num_machines` should be a number greater than 0.")

    @property
    def n_vcpus(self):
        return VCPUCount(
            self._cpu_info["cpu_cores_logical"] * self.num_machines,
            self._cpu_info["cpu_cores_logical"])

    def short_name(self) -> str:
        return "MachineGroup"

    @classmethod
    def from_api_response(cls, resp: dict):
        machine_group = super().from_api_response(resp)
        machine_group.num_machines = int(resp["max_vms"])
        machine_group.spot = bool(resp["spot"])
        return machine_group

    def __str__(self):
        return f"Machine Group {self.name} with {self.machine_type} machines"

    def _log_machine_group_info(self):
        super()._log_machine_group_info()
        logging.info("\t· Number of machines:         %s", self.num_machines)
        logging.info("\t· Spot:                       %s", self.spot)
        self.estimate_cloud_cost()


@dataclass(repr=False)
class ElasticMachineGroup(BaseMachineGroup):
    """Create an ElasticMachineGroup object.

    An ElasticMachineGroup is a set of identical machines that can
    automatically scale based on CPU load. The group starts with a
    minimum number of machines and adjusts its size as needed scaling
    to the maximum number of machines, ensuring both optimal performance
    and cost efficiency.
    Note: The machine group becomes active after calling the 'start' method,
    and billing commences once the machines are initiated.

    Args:
        machine_type: The type of GC machine to launch. Ex: "e2-standard-4".
            Check https://cloud.google.com/compute/docs/machine-resource for
        more information about machine types.
        zone: The zone where the machines will be launched.
        provider: The cloud provider of the machine group.
        threads_per_core: The number of threads per core (1 or 2).
        data_disk_gb: The size of the disk for user data (in GB).
        max_idle_time: Time without executing any task, after which the
            resource will be terminated.
        auto_terminate_ts: Moment in which the resource will be
            automatically terminated.
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
        auto_terminate_minutes: Duration, in minutes, the MPICluster will be
                kept alive. After auto_terminate_minutes minutes the machine
                will be terminated. This time will start counting after calling
                this method.
        min_machines: The minimum number of available machines. This is
            a quantity of machines that will be started initially and the
            minimum available machines, even in cases of low CPU load.
        max_machines: The maximum number of machines a machine group
            can scale up to.
        spot: Whether to use spot machines.
    """
    # Constructor arguments
    auto_resize_disk_max_gb: Optional[int] = None
    min_machines: int = 1
    max_machines: int = 2
    spot: bool = True

    # Internal attributes
    _is_elastic = True

    def __post_init__(self):
        """Validate inputs and initialize additional attributes after
        dataclass initialization."""
        super().__post_init__()

        self._active_machines = self.min_machines

        self._register_machine_group(min_vms=self.min_machines,
                                     max_vms=self.max_machines,
                                     is_elastic=self._is_elastic,
                                     num_vms=self._active_machines,
                                     spot=self.spot)

    def _validate_inputs(self):
        super()._validate_inputs()
        if self.min_machines < 0:
            raise ValueError(
                "`min_machines` should be a number equal or greater than 0.")

        if self.min_machines >= self.max_machines:
            raise ValueError("`max_machines` should be greater "
                             "than `min_machines`.")

    @property
    def n_vcpus(self):
        return VCPUCount(
            self._cpu_info["cpu_cores_logical"] * self.max_machines,
            self._cpu_info["cpu_cores_logical"])

    def short_name(self) -> str:
        return "ElasticMachineGroup"

    @classmethod
    def from_api_response(cls, resp: dict):
        machine_group = super().from_api_response(resp)
        machine_group.spot = bool(resp["spot"])
        machine_group.max_machines = int(resp["max_vms"])
        machine_group.min_machines = int(resp["min_vms"])
        return machine_group

    def active_machines_to_str(self) -> str:
        """Returns a string representation of the
        number of machines currently running.
        """
        return f"{self._active_machines}/{self.max_machines} (max)"

    def __str__(self):
        return f"Elastic Machine Group {self.name} with {self.machine_type} " \
             "machines"

    def _log_machine_group_info(self):
        super()._log_machine_group_info()
        logging.info("\t· Maximum number of machines: %s", self.max_machines)
        logging.info("\t· Minimum number of machines: %s", self.min_machines)
        logging.info("\t· Spot:                       %s", self.spot)
        self.estimate_cloud_cost()


@dataclass(repr=False)
class MPICluster(BaseMachineGroup):
    """Create a MPICluster object.

   A MPI cluster is a collection of homogenous machines all working together on
    a common task given the configurations that are launched in Google Cloud.
    Note: The cluster will be available only after calling 'start' method.
    The billing will start only after the machines are started.

    Args:
        machine_type: The type of GC machine to launch. Ex: "e2-standard-4".
            Check https://cloud.google.com/compute/docs/machine-resource for
            information about machine types.
        zone: The zone where the machines will be launched.
        provider: The cloud provider of the machine group.
        threads_per_core: The number of threads per core (1 or 2).
        data_disk_gb: The size of the disk for user data (in GB).
        max_idle_time: Time without executing any task, after which the
            resource will be terminated.
        auto_terminate_minutes: Duration, in minutes, the MPICluster will be
                kept alive. After auto_terminate_minutes minutes the machine
                will be terminated. This time will start counting after calling
                this method.
        auto_terminate_ts: Moment in which the resource will be
            automatically terminated.
        num_machines: The number of virtual machines to launch.
    """
    # Constructor arguments
    num_machines: int = 2
    spot: bool = True

    # Internal attributes
    auto_resize_disk_max_gb = None
    _type = ResourceType.MPI.value
    _is_elastic = False

    def __post_init__(self):
        """Validate inputs and initialize additional attributes after
        dataclass initialization."""
        super().__post_init__()

        self._register_machine_group(num_vms=self.num_machines,
                                     is_elastic=self._is_elastic,
                                     spot=self.spot,
                                     type=self._type)

    def _validate_inputs(self):
        super()._validate_inputs()
        if self.num_machines < 1:
            raise ValueError(
                "`num_machines` should be a number greater than 0.")

    @property
    def n_vcpus(self):
        return VCPUCount(
            self._cpu_info["cpu_cores_logical"] * self.num_machines,
            self._cpu_info["cpu_cores_logical"])

    @property
    def available_vcpus(self):
        """Returns the number of vCPUs available to the resource.

        For a mpi cluster with 2 machines, each with 4 vCPUs, this will
        return 8.
        """

        return self.n_vcpus.total

    def short_name(self) -> str:
        return "MPICluster"

    @classmethod
    def from_api_response(cls, resp: dict):
        machine_group = super().from_api_response(resp)
        machine_group.num_machines = int(resp["max_vms"])
        return machine_group

    def __str__(self):
        return f"MPI Cluster {self.name} with {self.machine_type} " \
               f"x{self.num_machines} machines"

    def _log_machine_group_info(self):
        super()._log_machine_group_info()
        logging.info("\t· Number of machines:       %s", self.num_machines)
        logging.info("\t· Spot:                     %s", self.spot)
        self.estimate_cloud_cost()

    def _log_estimated_spot_vm_savings(self) -> None:
        return


def _fetch_machine_groups_from_api():
    """Get all active machine groups of a user from the API."""
    try:
        api_compute = compute_api.ComputeApi(inductiva.api.get_client())
        response = api_compute.list_active_user_instance_groups()

        return json.loads(response.response.data)

    except inductiva.client.ApiException as api_exception:
        raise api_exception


def _get_machine_group_class(machine_type: str, is_elastic: bool):
    """Returns the class of the machine group"""
    if is_elastic:
        mg_class = ElasticMachineGroup
    elif machine_type == "standard":
        mg_class = MachineGroup
    elif machine_type == "mpi":
        mg_class = MPICluster
    else:
        raise ValueError("Unknown resource configuration.")
    return mg_class


def get_by_name(machine_name: str):
    """Returns the machine group corresponding to `machine_name`."""
    try:
        api_compute = compute_api.ComputeApi(inductiva.api.get_client())
        response = api_compute.get_vm_group_by_name({"name": machine_name})
        response = json.loads(response.response.data)
        mg_class = _get_machine_group_class(response["type"],
                                            response["is_elastic"])
        return mg_class.from_api_response(response)
    except inductiva.client.ApiException as api_exception:
        raise api_exception


def get():
    """Returns a list of 'Resource' objects."""

    # Retrive the active resource names
    machine_groups = _fetch_machine_groups_from_api()
    machine_group_list = []

    for mg in machine_groups:
        mg_class = _get_machine_group_class(mg["type"], mg["is_elastic"])
        machine_group_list.append(mg_class.from_api_response(mg))

    return machine_group_list
