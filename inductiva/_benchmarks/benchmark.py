"""API for Benchmarking"""
import enum
import json
import csv
from typing import Optional, Self, Union
from inductiva import types
from inductiva.simulators.simulator import Simulator
from inductiva.projects.project import Project


class ExportFormat(enum.Enum):
    """Enumeration of supported benchmark export formats."""
    JSON = "json"
    CSV = "csv"


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

    def export(
        self,
        fmt: Union[ExportFormat, str] = ExportFormat.JSON,
        filename: Optional[str] = None,
    ):
        """
        Exports the benchmark performance metrics in the specified format.

        Args:
            fmt (ExportFormat): The format to export the results in. Defaults
                to ExportFormat.JSON.
            filename (Optional[str]): The name of the output file to save the
                exported results. Defaults to the benchmark's name if not
                provided.
        """
        if isinstance(fmt, str):
            fmt = ExportFormat[fmt.upper()]
        metrics = self.gather_metrics()
        filename = filename or f"{self.name}.{fmt.value}"
        if fmt == ExportFormat.JSON:
            with open(filename, mode="w", encoding="utf-8") as file:
                json_content = json.dumps(obj=metrics, indent=4)
                file.write(json_content)
        elif fmt == ExportFormat.CSV:
            with open(filename, mode="w", encoding="utf-8") as file:
                csv_content = [{
                    "machine_type": machine_type,
                    **task,
                } for machine_type, tasks in metrics.items() for task in tasks]
                fieldnames = csv_content[0].keys()
                writer = csv.DictWriter(file, fieldnames=fieldnames)
                writer.writeheader()
                writer.writerows(csv_content)
        else:
            raise ValueError(f"Unsupported export format: {fmt}")

    def gather_metrics(self) -> dict:
        """
        Gathers performance metrics for all tasks associated with the
        benchmark, which include computation cost and execution time.    
    
        Returns:
            dict: A dictionary organized by virtual machine type, containing
                performance metrics for each task associated with the benchmark.
        """
        data = {}
        tasks = self.get_tasks()
        for task in tasks:
            info = task.get_info()
            vm_type = info.executer.vm_type
            data.setdefault(vm_type, [])
            data[vm_type].append({
                "task_id": info.task_id,
                "estimated_computation_cost": info.estimated_computation_cost,
                "computation_time": info.time_metrics.computation_seconds.value,
            })
        return data

    def gather_detailed_data(self) -> dict:
        """
        Gathers comprehensive information about all tasks associated with the 
        benchmark, including additional performance metrics and task metadata.

        Returns:
            dict: A dictionary organized by virtual machine type, containing
                detailed information about each task associated with the
                benchmark.
        """
        data = {}
        tasks = self.get_tasks()
        for task in tasks:
            info = task.get_info()
            vm_type = info.executer.vm_type
            data.setdefault(vm_type, [])
            input_filename = "input.json"
            input_dir_path = task.download_inputs(filenames=[input_filename])
            input_file_path = input_dir_path.joinpath(input_filename)
            with open(input_file_path, mode="r", encoding="utf-8") as file:
                input_json = json.load(file)
            data[vm_type].append({
                "task_info": info.to_dict(),
                "task_input_metadata": input_json,
            })
        return data
