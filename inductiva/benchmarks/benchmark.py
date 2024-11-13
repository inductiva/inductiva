"""API for Benchmarking"""
import enum
import json
import csv
from typing import Optional, Union
from typing_extensions import Self
from inductiva import types, resources
from inductiva.simulators import Simulator
from inductiva.projects import Project
from collections import defaultdict


class ExportFormat(enum.Enum):
    """Enumeration of supported benchmark export formats."""
    JSON = "json"
    CSV = "csv"


class ColumnExportMode(enum.Enum):
    """Enumeration of benchmark column export modes."""
    ALL = "all"
    DISTINCT = "distinct"


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
        self.kwargs = {**self.kwargs, **kwargs}
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
                machine_group.start(wait_on_pending_quota=False)
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
        columns: Union[ColumnExportMode, str] = ColumnExportMode.DISTINCT,
    ):
        """
        Exports the benchmark performance metrics in the specified format.

        Args:
            fmt (Union[ExportFormat, str]): The format to export the results
                in. Defaults to ExportFormat.JSON.
            filename (Optional[str]): The name of the output file to save the
                exported results. Defaults to the benchmark's name if not
                provided.
            columns (Union[ColumnExportMode, str]): The columns to include in
                the exported results. Defaults to ColumnExportMode.DISTINCT that
                includes only the parameters that vary between different runs.
        """
        if isinstance(fmt, str):
            fmt = ExportFormat[fmt.upper()]
        info = self.runs_info(columns=columns)
        filename = filename or f"{self.name}.{fmt.value}"
        if fmt == ExportFormat.JSON:
            with open(filename, mode="w", encoding="utf-8") as file:
                json_content = json.dumps(obj=info, indent=4)
                file.write(json_content)
        elif fmt == ExportFormat.CSV:
            with open(filename, mode="w", encoding="utf-8") as file:
                fieldnames = info[0].keys() if info else []
                writer = csv.DictWriter(file, fieldnames=fieldnames)
                writer.writeheader()
                writer.writerows(info)
        else:
            raise ValueError(f"Unsupported export format: {fmt}")

    def runs_info(
        self,
        columns: Union[ColumnExportMode, str] = ColumnExportMode.DISTINCT,
    ) -> list:
        """
        Gathers the configuration and performance metrics for each run
        associated with the benchmark in a list, including computation cost and
        execution time.
        
        Args:
            columns (Union[ColumnExportMode, str]): The columns to include in
                the exported results. Defaults to ColumnExportMode.DISTINCT that
                includes only the parameters that vary between different runs.

        Returns:
            list: A list containing the configuration and performance 
                metrics for each run.
        """

        def get_task_input_params(task):
            input_filename = "input.json"
            input_dir_path = task.download_inputs(filenames=[input_filename])
            input_file_path = input_dir_path.joinpath(input_filename)
            with open(input_file_path, mode="r", encoding="utf-8") as file:
                return json.load(file)

        def filter_distinct_columns(info):
            attrs_lsts = defaultdict(list)
            for attrs in info:
                for attr, value in attrs.items():
                    attrs_lsts[attr].append(value)
            filtered = {attr for attr, values in attrs_lsts.items() \
                        if values and len(values) != values.count(values[0])}
            return [{attr: attrs[attr] for attr in filtered} for attrs in info]

        if isinstance(columns, str):
            columns = ColumnExportMode[columns.upper()]

        info = []
        tasks = self.get_tasks()
        for task in tasks:
            task_input_params = get_task_input_params(task)
            task_info = task.info
            task_machine_type = task_info.executer.vm_type \
                if task_info.executer else None
            task_time = task_info.time_metrics.computation_seconds.value
            task_cost = task_info.estimated_computation_cost
            info.append({
                "task_id": task_info.task_id,
                "simulator": task_info.simulator,
                "machine_type": task_machine_type,
                "computation_time": task_time,
                "estimated_computation_cost": task_cost,
                **task_input_params,
            })

        return filter_distinct_columns(info) \
            if columns == ColumnExportMode.DISTINCT \
            else info

    def terminate(self) -> Self:
        """
        Terminates all active machine groups associated with the benchmark.
        """
        tasks = self.get_tasks()
        machines = {
            task.info.executer.uuid:
            resources.machine_groups.get_by_name(task.info.executer.vm_name)
            for task in tasks if task.info.executer
        }
        for machine in machines.values():
            machine.terminate(verbose=False)