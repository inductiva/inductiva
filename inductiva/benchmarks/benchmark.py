import datetime
import time
from typing import Any, Dict, List

import inductiva
from inductiva import types
from inductiva.benchmarks.benchmark_logger import BenchmarkLogger
from inductiva.benchmarks.to_json import ToJson


class Benchmark:
    """Class to create a benchmark."""

    def __init__(self,
                 name: str,
                 input_files: types.Path,
                 replicas: int,
                 machines: List[Dict[str, Any]],
                 simulator: inductiva.simulators.Simulator,
                 logger: BenchmarkLogger = ToJson(path="metadata"),
                 append: bool = False,
                 input_args: Dict[str, Any] = None,
                 machine_args: Dict[str, Any] = None,
                 simulator_args: Dict[str, Any] = None):
        self.name = name
        self.input_files = input_files
        self.replicas = replicas
        self.machines = machines
        self.simulator = simulator
        self.logger = logger
        self.logger.create_json(name)
        self.append = append
        self.input_args = input_args
        self.machine_args = machine_args
        self.simulator_args = simulator_args

        total_vcpu, total_cost, total_machines, total_tasks = self.report()

        self.total_vcpu = total_vcpu
        self.total_cost = total_cost
        self.total_machines = total_machines
        self.total_tasks = total_tasks

    def start(self):
        """Start the benchmark."""
        #TODO

        # Add docstring
        # write to csv
        # Test
        # Test
        # And more tests

        project = inductiva.projects.Project(self.name, append=True)
        project.open()

        tasks = project.get_tasks()
        if len(tasks) > 0 and not self.append:
            print(f"\nBenchmark {self.name} already exists.")
            #prompt the user if he wants to append the tasks to the existing benchmark
            #or cancel the operation
            print("Do you want to append the tasks to the existing benchmark?")
            print(
                "Type 'yes' to append the tasks or 'no' to cancel the operation."
            )
            print("You can also use the append argument to append the tasks to "
                  "the existing benchmark without any input necessary.")
            response = input()
            if response in ("no", "n"):
                print("Operation canceled.")
                project.close()
                return
            elif response in ("yes", "y"):
                print("Appending tasks to the existing benchmark.")

        quotas = inductiva.users.get_quotas()

        if (quotas["total_num_vcpus"]["in_use"] + self.total_vcpu
                > quotas["total_num_vcpus"]["max_allowed"] or
                quotas["total_num_machines"]["in_use"] + self.total_machines
                > quotas["total_num_machines"]["max_allowed"] or
                quotas["cost_per_hour"]["in_use"] + self.total_cost
                > quotas["cost_per_hour"]["max_allowed"]):
            print("This benchmark will exceed the current quotas.")
            print("Do you want to continue? (yes/no)")
            response = input()
            response = response.lower()

            if response == "no" or response == "n":
                print("Operation canceled.")
                return
            elif response == "yes" or response == "y":
                print("Continuing with the benchmark.")

        for machine in self.machines:

            #Iterate over machine_args and check if any value is callable,
            #if it is call it with machine["machine_type"] as parameter
            if self.machine_args is not None:
                for key, value in self.machine_args.items():
                    if callable(value):
                        self.machine_args[key] = value(machine["machine_type"])

            print(
                f"\nCreating machine group {machine['machine_type']} with {self.machine_args}"
            )

            machine_group = inductiva.resources.MachineGroup(
                machine_type=machine["machine_type"], **self.machine_args)

            #Iterate over input_args and simulator_args check if any value is callable,
            #if it is call it with machine_group as parameter
            if self.simulator_args is not None:
                for key, value in self.simulator_args.items():
                    if callable(value):
                        self.simulator_args[key] = value(machine_group)
            if self.input_args is not None:
                for key, value in self.input_args.items():
                    if callable(value):
                        self.input_args[key] = value(machine_group)

                template_manager = inductiva.TemplateManager(
                    template_dir=self.input_files)

                template_manager.set_root_dir(self.name)
                input_files = template_manager.render_dir(**self.input_args)
            else:
                input_files = self.input_files

            # input_files either come from the template manager or are passed directly
            self.simulator_args["input_dir"] = input_files
            self.simulator_args["on"] = machine_group

            print(
                f"\nRunning {self.simulator.api_method_name} with {self.simulator_args}"
            )
            if machine_group is not None:
                while not self._can_start_resource(machine_group):
                    print("This machine will exceed the current quotas.\n"
                          "Going to sleep and trying again later.")
                    time.sleep(60)

                machine_group.start(max_idle_time=datetime.timedelta(minutes=1))

            for _ in range(self.replicas):
                self.simulator.run(**self.simulator_args)

        project.close()

    def _can_start_resource(
            self, resource: inductiva.resources.machines_base.BaseMachineGroup
    ) -> bool:
        """Check if the resource can be started.
        
        This method checks if the resource can be started by checking
        the available quotas and resource usage.

        returns:
            bool: True if the resource can be started, False otherwise.
        """
        quotas = inductiva.users.get_quotas()
        cost_in_use = quotas["cost_per_hour"]["in_use"]
        vcpu_in_use = quotas["total_num_vcpus"]["in_use"]
        machines_in_use = quotas["total_num_machines"]["in_use"]

        cost_max = quotas["cost_per_hour"]["max_allowed"]
        vcpu_max = quotas["total_num_vcpus"]["max_allowed"]
        machines_max = quotas["total_num_machines"]["max_allowed"]

        current_vcpu = int(resource.machine_type.split("-")[2])

        if (cost_in_use + resource.estimate_cloud_cost() > cost_max or
                vcpu_in_use + current_vcpu > vcpu_max or
                machines_in_use + 1 > machines_max):
            return False
        return True

    def report(self, verbose: bool = True):
        """Prints a report of the benchmark.
        This includes:
            Total number of vcpus
            Total cost
            Total number of machines
            Total number of tasks
        This information will be used to promp the user if a benchmark is
        created that will exceed the quotas. This function will print the 
        information and return the values.
        Args:
            verbose (bool): If True, prints the report.
        return:
            Total number of vcpus
            Total cost
            Total number of machines
            Total number of tasks
        """
        total_vcpu = 0
        total_cost = 0.0
        total_machines = len(self.machines)
        total_number_tasks = self.replicas * total_machines

        for machine in self.machines:
            total_vcpu += machine["num_cpus"]
            total_cost += machine["price"]
            total_machines += 1
        if verbose:
            print(f"Report for benchmark: {self.name}")
            print(f"---> Total number of vcpus: {total_vcpu}")
            print(f"---> Total cost: {total_cost}")
            print(f"---> Total number of machines: {total_machines}")
            print(f"---> Total number of tasks: {total_number_tasks}")
        return total_vcpu, total_cost, total_machines, total_number_tasks
