"""Classes to manage different Google Cloud machine group types."""
from absl import logging

from inductiva.resources import machines_base


class MachineGroup(machines_base.BaseMachineGroup):
    """Class to launch and manage a group of machines in Google Cloud.

    A machine group is a collection of homogenous machines with given the
    configurations that are launched in Google Cloud.
    Note: The machine group will be available only after calling 'start' method.
    The billing will start only after the machines are started."""

    def __init__(
        self,
        machine_type: str,
        num_machines: int = 1,
        spot: bool = False,
        data_disk_gb: int = 10,
        register: bool = True,
    ) -> None:
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
            num_machines: The number of virtual machines to launch.
            spot: Whether to use spot machines.
            data_disk_gb: The size of the disk for user data (in GB).
        """
        super().__init__(machine_type=machine_type,
                         data_disk_gb=data_disk_gb,
                         register=register)
        self.num_machines = num_machines
        self.spot = spot
        self.__is_elastic = False

        if register:
            logging.info("Registering MachineGroup configurations:")
            self._register_machine_group(num_vms=self.num_machines,
                                         spot=self.spot,
                                         is_elastic=self.__is_elastic)

    @classmethod
    def from_api_response(cls, resp: dict):
        machine_group = super().from_api_response(resp)
        machine_group.num_machines = int(resp["num_vms"])
        machine_group.spot = bool(resp["spot"])
        machine_group.register = False
        return machine_group

    def __repr__(self):
        class_name = self.__class__.__name__
        return f"{class_name}(name=\"{self.name}\")"

    def __str__(self):
        return f"Machine Group {self.name} with {self.machine_type} machines"

    def start(self):
        """Starts all machines of the machine group."""
        return super().start(num_vms=self.num_machines,
                             is_elastic=self.__is_elastic,
                             spot=self.spot)

    def terminate(self):
        """Terminates all machines of the machine group."""
        return super().terminate(num_vms=self.num_machines,
                                 is_elastic=self.__is_elastic,
                                 spot=self.spot)

    def _log_machine_group_info(self):
        super()._log_machine_group_info()
        logging.info("> Number of machines: %s", self.num_machines)
        logging.info("> Spot:               %s", self.spot)
        self.estimate_cloud_cost()

    def estimate_cloud_cost(self):
        """Estimates a cost per hour of the machine group in US dollars.

        This is only an estimate of having a machine group with the
        specified configurations up in the cloud. The actual cost may vary.

        Returns:
            The estimated cost per hour of the machine group, in US
              dollars ($/h)."""
        cost_per_machine = super()._get_estimated_cost(self.spot)
        estimated_cost = cost_per_machine * self.num_machines
        logging.info("> Estimated cloud cost of machine group: %.3f $/h",
                     estimated_cost)
        return estimated_cost


class ElasticMachineGroup(machines_base.BaseMachineGroup):
    """Manages an elastic machine group in Google Cloud.

    An ElasticMachineGroup is a set of identical machines that can
    automatically scale based on CPU load. The group starts with a
    minimum number of machines and adjusts its size as needed scaling
    to the maximum number of machines, ensuring both optimal performance
    and cost efficiency.

    Note: The machine group becomes active after calling the 'start' method,
    and billing commences once the machines are initiated.
    """

    def __init__(
        self,
        machine_type: str,
        min_machines: int = 1,
        max_machines: int = 1,
        spot: bool = False,
        data_disk_gb: int = 10,
        register: bool = True,
    ) -> None:
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
            min_machines: The minimum number of available machines. This is
              a qunatity of machines that will be started initially and the
              minimum available machines, even in cases of low CPU load.
            max_machines: The maximum number of machines a machine group
              can scale up to.
            spot: Whether to use spot machines.
            data_disk_gb: The size of the disk for user data (in GB).
        """
        if max_machines < min_machines:
            raise ValueError("`max_machines` should be greater "
                             "than `min_machines`.")
        super().__init__(machine_type=machine_type,
                         data_disk_gb=data_disk_gb,
                         register=register)
        self.min_machines = min_machines
        self.max_machines = max_machines
        self.num_active_machines = min_machines
        self.__is_elastic = True
        self.spot = spot

        if self.register:
            logging.info("Registering ElasticMachineGroup configurations:")
            self._register_machine_group(min_vms=self.min_machines,
                                         max_vms=self.max_machines,
                                         is_elastic=self.__is_elastic,
                                         num_vms=self.num_active_machines,
                                         spot=self.spot)

    @classmethod
    def from_api_response(cls, resp: dict):
        machine_group = super().from_api_response(resp)
        machine_group.spot = bool(resp["spot"])
        machine_group.max_machines = int(resp["max_vms"])
        machine_group.min_machines = int(resp["min_vms"])
        machine_group.num_active_machines = int(resp["num_vms"])
        return machine_group

    def __repr__(self):
        class_name = self.__class__.__name__
        return f"{class_name}(name=\"{self.name}\")"

    def __str__(self):
        return f"Elastic Machine Group {self.name} with {self.machine_type} " \
             "machines"

    def start(self):
        """Starts minimum number of machines."""
        return super().start(num_vms=self.min_machines,
                             min_vms=self.min_machines,
                             max_vms=self.max_machines,
                             is_elastic=self.__is_elastic,
                             spot=self.spot)

    def terminate(self):
        """Terminates all machines of the machine group."""
        return super().terminate(num_vms=self.min_machines,
                                 min_vms=self.min_machines,
                                 max_vms=self.max_machines,
                                 is_elastic=self.__is_elastic,
                                 spot=self.spot)

    def _log_machine_group_info(self):
        super()._log_machine_group_info()
        logging.info("> Maximum number of machines: %s", self.max_machines)
        logging.info("> Minimum number of machines: %s", self.min_machines)
        logging.info("> Spot: %s", self.spot)
        self.estimate_cloud_cost()

    def estimate_cloud_cost(self):
        """Estimates a cost per hour of min and max machines in US dollars.

        these are the estimted costs of having minimum and the
        maximum number of machines up in the cloud. The final cost will vary
        depending on the total usage of the machines."""
        cost_per_machine = super()._get_estimated_cost(self.spot)
        logging.info(
            "Note: these are the estimated costs of having minimum and the "
            "maximum number of machines up in the cloud. The final cost will "
            "vary depending on the total usage of the machines.")
        logging.info(
            "> Minimum estimated cloud cost of elastic machine group: "
            "%.3f $/h.", cost_per_machine * self.min_machines)
        logging.info(
            "> Maximum estimated cloud cost  of elastic machine group:"
            " %.3f $/h.", cost_per_machine * self.max_machines)
