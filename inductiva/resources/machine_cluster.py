"""Class to manage the MPI cluster in Google Cloud."""
from absl import logging

from inductiva.resources import machines_base


class MPICluster(machines_base.BaseMachineGroup):
    """Class to launch and manage an MPI cluster in Google Cloud.

    A MPI cluster is a collection of homogenous machines all working together on
    a common task given the configurations that are launched in Google Cloud.
    Note: The cluster will be available only after calling 'start' method.
    The billing will start only after the machines are started."""

    def __init__(self,
                 machine_type: str,
                 num_machines: int = 2,
                 data_disk_gb: int = 10,
                 register: bool = True) -> None:
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
            data_disk_gb: The size of the disk for user data (in GB).
        """
        super().__init__(machine_type=machine_type,
                         data_disk_gb=data_disk_gb,
                         register=register)
        self.num_machines = num_machines
        self.__type = machines_base.ResourceType.MPI.value
        self.__is_elastic = False
        self.__spot = False

        if register:
            logging.info("Registering MPICluster configurations:")
            self._register_machine_group(num_vms=self.num_machines,
                                         is_elastic=self.__is_elastic,
                                         spot=self.__spot,
                                         type=self.__type)

    @classmethod
    def from_api_response(cls, resp: dict):
        machine_group = super().from_api_response(resp)
        machine_group.num_machines = int(resp["num_vms"])
        machine_group.register = False
        return machine_group

    def __repr__(self):
        class_name = self.__class__.__name__
        return f"{class_name}(name=\"{self.name}\")"

    def __str__(self):
        return f"MPI Cluster {self.name} with {self.machine_type} " \
               f"x{self.num_machines} machines"

    def start(self):
        """Start the MPI Cluster."""
        return super().start(num_vms=self.num_machines,
                             is_elastic=self.__is_elastic,
                             spot=self.__spot,
                             type=self.__type)

    def terminate(self):
        """Terminates the MPI Cluster."""
        return super().terminate(num_vms=self.num_machines,
                                 is_elastic=self.__is_elastic,
                                 spot=self.__spot)

    def _log_machine_group_info(self):
        super()._log_machine_group_info()
        logging.info("> Number of machines: %s", self.num_machines)
        self.estimate_cloud_cost()

    def estimate_cloud_cost(self):
        """Estimates a cost per hour of the MPI cluster in US dollars.

        This is an estimate of the cost of MPI cluster with the
        specified configurations up in the cloud. The actual cost may vary.

        Returns:
            The estimated cost per hour of the machine group, in US
              dollars ($/h)."""
        #TODO: Contemplate disk size in the price.
        estimated_cost = super()._get_estimated_cost() * self.num_machines
        logging.info("> Estimated cloud cost of the MPI cluster: %.3f $/h",
                     estimated_cost)
        return estimated_cost
