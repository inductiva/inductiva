"""API for Benchmarking"""
import enum
import json
import csv
import logging
from typing import Optional, Union
from typing_extensions import Self
from collections import defaultdict
from inductiva import types, resources, projects, simulators, client
from inductiva.client.models import TaskStatusCode
from inductiva.utils.format_utils import CURRENCY_SYMBOL, TIME_UNIT


class ExportFormat(enum.Enum):
    """Enumeration of supported benchmark export formats."""
    JSON = "json"
    CSV = "csv"


class SelectMode(enum.Enum):
    """
    Enumeration of supported data selection modes, specifying which data 
    should be included in the benchmarking results.
    """
    ALL = "all"
    DISTINCT = "distinct"


class Benchmark(projects.Project):
    """Represents the benchmark runner."""

    class InfoKey:
        """Constant keys for the info used in the benchmark runs.

        :meta private:
        """
        TASK_ID = "task_id"
        SIMULATOR = "simulator"
        MACHINE_TYPE = "machine_type"
        TIME = f"computation_time ({TIME_UNIT})"
        COST = f"estimated_computation_cost ({CURRENCY_SYMBOL})"

    def __init__(self, name: str):
        """
        Initializes a new Benchmark instance.

        Args:
            name (str): The name of the benchmark runner.
        """
        super().__init__(name=name)
        self.runs = []
        self.simulator = None
        self.input_dir = None
        self.on = None
        self.kwargs = {}

    def set_default(
        self,
        simulator: Optional[simulators.Simulator] = None,
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
        simulator: Optional[simulators.Simulator] = None,
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

    def run(self, num_repeats: int = 2, wait_for_quotas: bool = True) -> Self:
        """
        Executes all added runs.

        Each run is executed the specified number of times, and the collection
        of runs is cleared afterwards.

        Args:
            num_repeats (int): The number of times to repeat each simulation
                run (default is 2).
            wait_for_quotas (bool): Indicates whether to wait for quotas to 
                become available before starting each resource. If `True`, the 
                program will actively wait in a loop, periodically sleeping and 
                checking for quotas. If `False`, the program crashes if quotas 
                are not available (default is `True`).

        Returns:
            Self: The current instance for method chaining.
        """
        for simulator, input_dir, machine_group, kwargs in self.runs:
            if not machine_group.started:
                machine_group.start(wait_for_quotas=wait_for_quotas)
            for _ in range(num_repeats):
                simulator.run(input_dir=input_dir,
                              on=machine_group,
                              project=self.name,
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
        status: Optional[Union[TaskStatusCode, str]] = None,
        select: Union[SelectMode, str] = SelectMode.DISTINCT,
    ):
        """
        Exports the benchmark performance metrics in the specified format.

        Args:
            fmt (Union[ExportFormat, str]): The format to export the results
                in. Defaults to ExportFormat.JSON.
            filename (Optional[str]): The name of the output file to save the
                exported results. Defaults to the benchmark's name if not
                provided.
            status (Optional[Union[TaskStatusCode, str]]): The status of the
                tasks to include in the benchmarking results. Defaults to None,
                which includes all tasks.
            select (Union[SelectMode, str]): The data to include in
                the benchmarking results. Defaults to SelectMode.DISTINCT that
                includes only the parameters that vary between different runs.
        """
        if isinstance(fmt, str):
            fmt = ExportFormat[fmt.upper()]
        info = self.runs_info(status=status, select=select)
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
        status: Optional[Union[TaskStatusCode, str]] = None,
        select: Union[SelectMode, str] = SelectMode.DISTINCT,
    ) -> list:
        """
        Gathers the configuration and performance metrics for each run
        associated with the benchmark in a list, including computation cost and
        execution time.
        
        Args:
            status (Optional[Union[TaskStatusCode, str]]): The status of the
                tasks to include in the benchmarking results. Defaults to None,
                which includes all tasks.
            select (Union[SelectMode, str]): The data to include in
                the benchmarking results. Defaults to SelectMode.DISTINCT that
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

        def select_distinct(info):
            attrs_lsts = defaultdict(list)
            for attrs in info:
                for attr, value in attrs.items():
                    attrs_lsts[attr].append(value)
            filtered = {attr for attr, values in attrs_lsts.items() \
                        if values and len(values) != values.count(values[0])}
            return [{attr: attrs[attr] for attr in filtered} for attrs in info]

        if isinstance(select, str):
            select = SelectMode[select.upper()]

        info = []
        tasks = self.get_tasks(status=status)
        for task in tasks:
            task_input_params = get_task_input_params(task)
            task_info = task.info
            task_machine_type = task_info.executer.vm_type \
                if task_info.executer else None
            task_time = task_info.time_metrics.computation_seconds.value
            task_cost = task_info.estimated_computation_cost
            info.append({
                Benchmark.InfoKey.TASK_ID: task_info.task_id,
                Benchmark.InfoKey.SIMULATOR: task_info.simulator,
                Benchmark.InfoKey.MACHINE_TYPE: task_machine_type,
                Benchmark.InfoKey.TIME: task_time,
                Benchmark.InfoKey.COST: task_cost,
                **task_input_params,
            })

        return select_distinct(info) \
            if select == SelectMode.DISTINCT \
            else info

    def terminate(self) -> Self:
        """
        Terminates all active machine groups associated with the
        benchmark.

        Returns:
            Self: The current instance for method chaining.
        """

        def _handle_suffix(executer):
            if executer.host_type == resources.utils.ProviderType.GCP:
                return "-".join(executer.vm_name.split("-")[:-1])
            return executer.vm_name

        machine_names = {
            _handle_suffix(task.info.executer)
            for task in self.get_tasks()
            if task.info.executer
        }

        for machine in resources.get():
            if machine.name not in machine_names:
                continue
            try:
                machine.terminate(verbose=False)
            except client.ApiException as api_exception:
                logging.warning(api_exception)

        return self
