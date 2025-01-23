"""Class to manage the MPI cluster in Google Cloud."""
import logging
from typing import Optional
import datetime

from inductiva.resources import machines_base


class MPICluster(machines_base.BaseMachineGroup):
    """Class to launch and manage an MPI cluster in Google Cloud.

    A MPI cluster is a collection of homogenous machines all working together on
    a common task given the configurations that are launched in Google Cloud.
    Note: The cluster will be available only after calling 'start' method.
    The billing will start only after the machines are started."""

    def __init__(
        self,
        machine_type: str,
        num_machines: int = 2,
        threads_per_core: int = 2,
        data_disk_gb: int = 10,
        max_idle_time: Optional[datetime.timedelta] = None,
        auto_terminate_ts: Optional[datetime.datetime] = None,
        register: bool = True,
    ) -> None:
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
            num_machines: The number of virtual machines to launch.
            threads_per_core: The number of threads per core (1 or 2).
            data_disk_gb: The size of the disk for user data (in GB).
            max_idle_time: Time without executing any task, after which the
              resource will be terminated.
            auto_terminate_ts: Moment in which the resource will be
              automatically terminated.
        """
        if num_machines < 1:
            raise ValueError(
                "`num_machines` should be a number greater than 0.")

        super().__init__(
            machine_type=machine_type,
            threads_per_core=threads_per_core,
            data_disk_gb=data_disk_gb,
            max_idle_time=max_idle_time,
            auto_terminate_ts=auto_terminate_ts,
            register=register,
        )

        # num_machines is the number of machines requested
        self.num_machines = num_machines
        #Number of active machines at the time of
        #the request machine_groups.get()
        self._active_machines = 0
        self.__type = machines_base.ResourceType.MPI.value
        self.__is_elastic = False
        self.__spot = False

        if register:
            logging.info("■ Registering MPICluster configurations:")
            self._register_machine_group(num_vms=self.num_machines,
                                         is_elastic=self.__is_elastic,
                                         spot=self.__spot,
                                         type=self.__type)

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

    def __repr__(self):
        class_name = self.__class__.__name__
        return f"{class_name}(name=\"{self.name}\")"

    def __str__(self):
        return f"MPI Cluster {self.name} with {self.machine_type} " \
               f"x{self.num_machines} machines"

    def start(self, wait_for_quotas: bool = False):
        """Start the MPI Cluster.
        Args:
            wait_for_quotas: If True, the method will wait for quotas to
              become available before starting the resource.
        """
        return super().start(
            wait_for_quotas=wait_for_quotas,
            is_elastic=self.__is_elastic,
            num_vms=self.num_machines,
            spot=self.__spot,
            type=self.__type,
        )

    def terminate(self, verbose: bool = True):
        """Terminates the MPI Cluster."""
        return super().terminate(num_vms=self.num_machines,
                                 is_elastic=self.__is_elastic,
                                 spot=self.__spot,
                                 verbose=verbose)

    def _log_machine_group_info(self):
        super()._log_machine_group_info()
        logging.info("\t· Number of machines:       %s", self.num_machines)
        self.estimate_cloud_cost()

    def _log_estimated_spot_vm_savings(self) -> None:
        return
