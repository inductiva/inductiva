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
        disk_size_gb: int = 70,
        zone: str = "europe-west1-b",
        register: bool = True,
    ) -> None:
        """Create a MachineGroup object.

        Args:
            machine_type: The type of GC machine to launch. Ex: "e2-standard-4".
              Check https://cloud.google.com/compute/docs/machine-resource for
              information about machine types.
            num_machines: The number of virtual machines to launch.
            spot: Whether to use spot machines.
            disk_size_gb: The size of the disk in GB, recommended min. is 60 GB.
            zone: The zone where the machines will be launched.
        """
        super().__init__(machine_type=machine_type,
                         spot=spot,
                         disk_size_gb=disk_size_gb,
                         zone=zone,
                         register=register)
        self.num_machines = num_machines
        self.is_elastic = False

        if register:
            super()._register_machine_group(num_instances=self.num_machines,
                                            is_elastic=self.is_elastic)
            self._log_machine_group_info()

    @classmethod
    def from_api_response(cls, resp: dict):
        machine_group = super().from_api_response(resp)
        machine_group.num_machines = resp["num_instances"]
        machine_group.register = False
        return machine_group

    def start(self):
        """Starts all machines of the machine group."""
        return super().start(num_instances=self.num_machines,
                             is_elastic=self.is_elastic)

    def terminate(self):
        """Terminates all machines of the machine group."""
        return super().terminate(num_instances=self.num_machines,
                                 is_elastic=self.is_elastic)

    def _log_machine_group_info(self):
        super()._log_machine_group_info()
        logging.info("> Number of machines: %s", self.num_machines)
        self.estimate_cloud_cost()

    def estimate_cloud_cost(self):
        """Estimates a cost per hour of the machine group in US dollars.

        This is only an estimate of having a machine group with the
        specified configurations up in the cloud. The actual cost may vary.

        Returns:
            The estimated cost per hour of the machine group in US
              dollars ($/h)."""
        #TODO: Contemplate disk size in the price.
        estimated_cost = round(
            float(super()._get_estimated_cost() * self.num_machines), 3)
        logging.info("Estimated cloud cost per hour for all machines : %s $/h",
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
        disk_size_gb: int = 70,
        zone: str = "europe-west1-b",
        register: bool = True,
    ) -> None:
        """Create an ElasticMachineGroup object.

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
            disk_size_gb: The size of the disk in GB, recommended min. is 60 GB.
            zone: The zone where the machines will be launched.
        """
        if max_machines < min_machines:
            raise ValueError("`max_machines` should be greater "
                             "than `min_machines`.")
        super().__init__(machine_type=machine_type,
                         spot=spot,
                         disk_size_gb=disk_size_gb,
                         zone=zone,
                         register=register)
        self.min_machines = min_machines
        self.max_machines = max_machines
        self.num_active_machines = min_machines
        self.is_elastic = True

        if self.register:
            super()._register_machine_group(
                min_instances=self.min_machines,
                max_instances=self.max_machines,
                is_elastic=self.is_elastic,
                num_instances=self.num_active_machines)
            self._log_machine_group_info()

    @classmethod
    def from_api_response(cls, resp: dict):
        machine_group = super().from_api_response(resp)
        machine_group.max_machines = resp["max_instances"]
        machine_group.min_machines = resp["min_instances"]
        machine_group.num_active_machines = resp["num_instances"]
        return machine_group

    def start(self):
        """Starts minimum number of machines."""
        return super().start(num_instances=self.min_machines,
                             min_instances=self.min_machines,
                             max_instances=self.max_machines,
                             is_elastic=self.is_elastic)

    def terminate(self):
        """Terminates all machines of the machine group."""
        return super().terminate(num_instances=self.min_machines,
                                 min_instances=self.min_machines,
                                 max_instances=self.max_machines,
                                 is_elastic=self.is_elastic)

    def _log_machine_group_info(self):
        super()._log_machine_group_info()
        logging.info("> Maximum number of machines: %s", self.max_machines)
        logging.info("> Minimum number of machines: %s", self.min_machines)
        self.estimate_cloud_cost()

    def estimate_cloud_cost(self):
        """Estimates a cost per hour of min and max machines in US dollars.

        these are the estimted costs of having minimum and the
        maximum number of machines up in the cloud. The final cost will vary
        depending on the total usage of the machines."""
        cost = super()._get_estimated_cost()
        logging.info(
            "Note: these are the estimted costs of having minimum and the "
            "maximum number of machines up in the cloud. The final cost will "
            "vary depending on the total usage of the machines.")
        logging.info("Minimum estimated cloud cost: "
                     "%s $/h.", round(float(cost * self.min_machines), 3))
        logging.info("Maximum estimated cloud cost: "
                     "%s $/h.", round(float(cost * self.max_machines), 3))
