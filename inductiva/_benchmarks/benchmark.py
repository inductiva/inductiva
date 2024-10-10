import json
from typing import Self
from inductiva import types
from inductiva.simulators.simulator import Simulator
from inductiva.projects.project import Project
from inductiva.resources import machine_groups


class Benchmark(Project):
    """Represents the benchmark runner."""

    def __init__(self, name: str, append: bool = True):
        """
        Initializes a new Benchmark instance.

        Args:
            name (str): The name of the benchmark runner.
            append (bool): Indicates whether to allow adding runs to the 
            existing benchmark (default is True).
        """
        super().__init__(name=name, append=append)
        self.runs = []

    def add_run(
        self,
        simulator: Simulator,
        input_dir: str,
        on: types.ComputationalResources,
        **kwargs,
    ) -> Self:
        """
        Adds a simulation run to the benchmark.

        Args:
            simulator (Simulator): The simulator to be used for this run.
            input_dir (str): The directory containing input files for the
            simulation.
            on (types.ComputationalResources): The computational resources to
            run the simulation on.
            **kwargs: Additional keyword arguments for the simulator run.

        Returns:
            Self: The current instance for method chaining.
        """
        self.runs.append((simulator, input_dir, on, kwargs))
        return self

    def run(self, num_repeats: int = 2) -> Self:
        """
        Executes all added runs.

        Each run is executed the specified number of times, and the collection
        of runs is cleared afterwards.

        Args:
            num_repeats (int): The number of times to repeat each simulation
            run (default is 2).

        Returns:
            Self: The current instance for method chaining.
        """
        with self:
            for simulator, input_dir, machine_group, kwargs in self.runs:
                machine_group.start(wait_on_pending_quota=True)
                for _ in range(num_repeats):
                    simulator.run(input_dir=input_dir,
                                  on=machine_group,
                                  **kwargs)
            self.runs.clear()
        return self

    def wait(self) -> Self:
        """
        Waits for all running tasks to complete.

        Returns:
            Self: The current instance for method chaining.
        """
        tasks = self.get_tasks()
        for task in tasks:
            task.wait(download_std_on_completion=False)
        return self

    def to_dict(self):
        """
        Compiles results from all completed tasks into a dictionary.

        Returns:
            dict: A dictionary containing the results of each task, organized 
            by virtual machine type.
        """
        results = {}
        tasks = self.get_tasks()
        for task in tasks:
            task_info = task.info
            tmp1, tmp2, *_ = task_info.executer.vm_name.split("-")
            machine_group_name = tmp1 + "-" + tmp2
            resource = machine_groups.get_by_name(machine_group_name)
            vm_type = task_info.executer.vm_type
            if vm_type not in results:
                cost_per_hour = resource.estimate_cloud_cost(verbose=False)
                results[vm_type] = {
                    "cost per hour": cost_per_hour,
                    "machine type": resource.__class__.__name__,
                    "provider": task_info.executer.host_type,
                    "machine count": resource.num_machines,
                    "tasks": []
                }
            metadata_filename = "input.json"
            inputs_path = task.download_inputs(filenames=[metadata_filename])
            metadata_path = inputs_path.joinpath(metadata_filename)
            with open(metadata_path, mode="r", encoding="utf-8") as file:
                metadata_content = json.load(file)
            results[vm_type]["tasks"].append({
                "info": task_info.to_dict(),
                "metadata": metadata_content,
            })
        return results
