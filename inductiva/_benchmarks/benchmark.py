"""API for Benchmarking"""
import json
from typing import Optional, Self
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
        self.simulator = None
        self.input_dir = None
        self.on = None
        self.kwargs = {}


    def set_default(
        self,
        simulator: Optional[Simulator] = None,
        input_dir: Optional[str] = None,
        on: Optional[types.ComputationalResources] = None,
        **kwargs,
    ) -> Self:
        """
        Sets default parameters for the benchmark runner.

        This method allows you to configure default settings for the benchmark,
        which will be used in subsequent runs unless explicitly overridden.

        Args:
            simulator (Optional[Simulator]): The simulator instance to be used
                as the default for future runs. If not provided, the current
                simulator will remain unchanged.
            input_dir (Optional[str]): The directory path for input files. If
                not provided, the current input directory will remain unchanged.
            on (Optional[types.ComputationalResources]): The computational
                resources to use for running the simulations. If not specified,
                the current resources will remain unchanged.
            **kwargs: Additional keyword arguments to set as default parameters 
                for the simulations. These will update any existing parameters
                with the same names.

        Returns:
            Self: The current instance for method chaining.
        """
        self.simulator = simulator or self.simulator
        self.input_dir = input_dir or self.input_dir
        self.on = on or self.on
        self.kwargs = kwargs or self.kwargs
        return self

    def add_run(
        self,
        simulator: Optional[Simulator] = None,
        input_dir: Optional[str] = None,
        on: Optional[types.ComputationalResources] = None,
        **kwargs,
    ) -> Self:
        """
        Adds a simulation run to the benchmark.

        Args:
            simulator (Optional[Simulator]): The simulator to be used for this
                run. If not provided, the previously set simulator will be used.
            input_dir (Optional[str]): The directory containing input files for
                the simulation. If not provided, the previously set input
                directory will be used.
            on (Optional[types.ComputationalResources]): The computational
                resources to run the simulation on. If not provided, the
                previously set resources will be used.
            **kwargs: Additional keyword arguments for the simulator run. These
                will overwrite any previously set parameters.

        Returns:
            Self: The current instance for method chaining.
        """
        self.runs.append((
            simulator or self.simulator,
            input_dir or self.input_dir,
            on or self.on,
            {
                **self.kwargs,
                **kwargs
            },
        ))
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
