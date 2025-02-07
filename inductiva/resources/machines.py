"""Classes to manage different Google Cloud machine group types."""
from dataclasses import dataclass
import logging
from typing import Optional

from inductiva.resources import machines_base


@dataclass
class MachineGroup(machines_base.BaseMachineGroup):
    """Create a MachineGroup object.
    
    The register argument is used to indicate if the machine group should
    be registered or if it was already registered. If set as False on
    initialization, then, the machine group is not registered and it
    can not be started in the cloud. This serves has an helper argument for
    retrieving already registered machine groups that can be started, for
    example, when retrieving with the `machines_groups.get` method.
    Users should not set this argument.

    Args:
        machine_type: The type of GC machine to launch. Ex: "e2-standard-4".
          Check https://cloud.google.com/compute/docs/machine-resource for
          information about machine types.
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
        num_machines: The number of virtual machines to launch.
        spot: Whether to use spot machines.
    """
    """
    A machine group is a collection of homogenous machines with given the
    configurations that are launched in Google Cloud.
    
    Note: The machine group will be available only after calling 'start' method.
    The billing will start only after the machines are started.
    """
    # Construtor arguments
    auto_resize_disk_max_gb: Optional[int] = None
    num_machines: int = 1
    spot: bool = True

    # Internal arguments
    _is_elastic = False

    def __post_init__(self):
        """Validate inputs and initialize additional attributes after
        dataclass initialization."""
        if self.num_machines < 1:
            raise ValueError(
                "`num_machines` should be a number greater than 0.")

        super().__post_init__()

        if self.register:
            logging.info("■ Registering MachineGroup configurations:")
            self._register_machine_group(num_vms=self.num_machines,
                                         spot=self.spot,
                                         is_elastic=self._is_elastic)

    @classmethod
    def from_api_response(cls, resp: dict):
        machine_group = super().from_api_response(resp)
        machine_group.num_machines = int(resp["max_vms"])
        machine_group.provider = resp["provider_id"]
        machine_group.__dict__["_active_machines"] = int(resp["num_vms"])
        machine_group.__dict__["machines"] = resp["machines"]
        machine_group.spot = bool(resp["spot"])
        machine_group.register = False
        return machine_group

    def __str__(self):
        return f"Machine Group {self.name} with {self.machine_type} machines"

    def _log_machine_group_info(self):
        super()._log_machine_group_info()
        logging.info("\t· Number of machines:         %s", self.num_machines)
        logging.info("\t· Spot:                       %s", self.spot)
        self.estimate_cloud_cost()


@dataclass
class ElasticMachineGroup(machines_base.BaseMachineGroup):
    """Create an ElasticMachineGroup object.

    The register argument is used to indicate if the machine group should
    be registered or if it was already registered. If set as False on
    initialization, then, the machine group is not registered and it
    can not be started in the cloud. This serves has an helper argument for
    retrieving already registered machine groups that can be started, for
    example, when retrieving with the `machines_groups.get` method.
    Users should not set this argument.

    Args:
        machine_type: The type of GC machine to launch. Ex: "e2-standard-4".
            Check https://cloud.google.com/compute/docs/machine-resource for
        more information about machine types.
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
        min_machines: The minimum number of available machines. This is
            a quantity of machines that will be started initially and the
            minimum available machines, even in cases of low CPU load.
        max_machines: The maximum number of machines a machine group
            can scale up to.
        spot: Whether to use spot machines.
    """
    """
    An ElasticMachineGroup is a set of identical machines that can
    automatically scale based on CPU load. The group starts with a
    minimum number of machines and adjusts its size as needed scaling
    to the maximum number of machines, ensuring both optimal performance
    and cost efficiency.

    Note: The machine group becomes active after calling the 'start' method,
    and billing commences once the machines are initiated.
    """
    # Construtor arguments
    auto_resize_disk_max_gb: Optional[int] = None
    min_machines: int = 1
    max_machines: int = 2
    spot: bool = True

    # Internal arguments
    _is_elastic = True

    def __post_init__(self):
        """Validate inputs and initialize additional attributes after
        dataclass initialization."""
        if self.min_machines < 0:
            raise ValueError(
                "`min_machines` should be a number equal or greater than 0.")

        if self.min_machines >= self.max_machines:
            raise ValueError("`max_machines` should be greater "
                             "than `min_machines`.")

        super().__post_init__()

        self._active_machines = self.min_machines

        if self.register:
            logging.info("■ Registering ElasticMachineGroup configurations:")
            self._register_machine_group(min_vms=self.min_machines,
                                         max_vms=self.max_machines,
                                         is_elastic=self._is_elastic,
                                         num_vms=self._active_machines,
                                         spot=self.spot)

    @classmethod
    def from_api_response(cls, resp: dict):
        machine_group = super().from_api_response(resp)
        machine_group.spot = bool(resp["spot"])
        machine_group.max_machines = int(resp["max_vms"])
        machine_group.min_machines = int(resp["min_vms"])
        machine_group.__dict__["_active_machines"] = int(resp["num_vms"])
        machine_group.__dict__["machines"] = resp["machines"]
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


@dataclass
class MPICluster(machines_base.BaseMachineGroup):
    """Create a MPICluster object.

    The register argument is used to indicate if the machine group should
    be registered or if it was already registered. If set as False on
    initialization, then, the machine group is not registered and it
    can not be started in the cloud. This serves has an helper argument for
    retrieving already registered machine groups that can be started, for
    example, when retrieving with the `machines_groups.get` method.
    Users should not set this argument.

    Args:
        machine_type: The type of GC machine to launch. Ex: "e2-standard-4".
            Check https://cloud.google.com/compute/docs/machine-resource for
            information about machine types.
        provider: The cloud provider of the machine group.
        threads_per_core: The number of threads per core (1 or 2).
        data_disk_gb: The size of the disk for user data (in GB).
        max_idle_time: Time without executing any task, after which the
            resource will be terminated.
        auto_terminate_ts: Moment in which the resource will be
            automatically terminated.
        num_machines: The number of virtual machines to launch.
    """
    """
    A MPI cluster is a collection of homogenous machines all working together on
    a common task given the configurations that are launched in Google Cloud.
    
    Note: The cluster will be available only after calling 'start' method.
    The billing will start only after the machines are started.
    """
    # Construtor arguments
    num_machines: int = 2

    # Internal arguments
    auto_resize_disk_max_gb = None
    _type = machines_base.ResourceType.MPI.value
    _is_elastic = False
    _spot = False

    def __post_init__(self):
        """Validate inputs and initialize additional attributes after
        dataclass initialization."""
        if self.num_machines < 1:
            raise ValueError(
                "`num_machines` should be a number greater than 0.")

        super().__post_init__()

        if self.register:
            logging.info("■ Registering MPICluster configurations:")
            self._register_machine_group(num_vms=self.num_machines,
                                         is_elastic=self._is_elastic,
                                         spot=self._spot,
                                         type=self._type)

    @property
    def available_vcpus(self):
        """Returns the number of vCPUs available to the resource.

        For a mpi cluster with 2 machines, each with 4 vCPUs, this will
        return 8.
        """

        return self.n_vcpus.total

    @classmethod
    def from_api_response(cls, resp: dict):
        machine_group = super().from_api_response(resp)
        machine_group.num_machines = int(resp["max_vms"])
        machine_group.__dict__["_active_machines"] = int(resp["num_vms"])
        machine_group.__dict__["machines"] = resp["machines"]
        machine_group.register = False
        return machine_group

    def __str__(self):
        return f"MPI Cluster {self.name} with {self.machine_type} " \
               f"x{self.num_machines} machines"

    def _log_machine_group_info(self):
        super()._log_machine_group_info()
        logging.info("\t· Number of machines:       %s", self.num_machines)
        self.estimate_cloud_cost()

    def _log_estimated_spot_vm_savings(self) -> None:
        return
