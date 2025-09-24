"""API for Benchmarking"""
import enum
import json
import csv
import logging
from typing import Optional, Union
import tqdm
from tqdm.contrib.logging import logging_redirect_tqdm
from typing_extensions import Self
from inductiva import types, resources, projects, simulators, client
from inductiva.client.models import TaskStatusCode
from inductiva.projects.project import ProjectType
from inductiva.utils.format_utils import CURRENCY_SYMBOL, TIME_UNIT

from concurrent.futures import ThreadPoolExecutor, as_completed


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

    def __init__(self, name: str, verbose: bool = False):
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
        self.verbose = verbose

    def _get_project_type(self):
        return ProjectType.BENCHMARK

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
        if not self.verbose:
            logging.info(
                "Preparing tasks for Benchmark to start. "
                "This may take a few minutes.\n"
                "Note that stopping this process will \033[1minterrupt\033[0m "
                "the submission of the tasks. Please wait...\n")
        for simulator, input_dir, machine_group, kwargs in self.runs:
            if not machine_group.started:
                machine_group.start(wait_for_quotas=wait_for_quotas,
                                    verbose=self.verbose)
            for _ in range(num_repeats):
                simulator.run(input_dir=input_dir,
                              on=machine_group,
                              project=self.name,
                              resubmit_on_preemption=True,
                              verbose=self.verbose,
                              **kwargs)
        self.runs.clear()
        logging.info(
            "■ Benchmark \033[1m%s\033[0m has started.\n"
            "  Go to https://console.inductiva.ai/benchmarks/%s "
            "for more details.\n", self.name, self.name)
        return self

    def wait(self) -> Self:
        """
        Waits for all running tasks to complete.

        Returns:
            Self: The current instance for method chaining.
        """
        tasks = self.get_tasks()

        completed_tasks = [task for task in tasks if task.info.is_terminal]
        running_tasks = [task for task in tasks if not task.info.is_terminal]

        logging.info("Waiting for Benchmark \033[1m%s\033[0m to complete...\n",
                     self.name)

        with tqdm.tqdm(total=len(tasks),
                       desc="Running Benchmark",
                       bar_format="{l_bar}{bar}| {n_fmt}/{total_fmt} {unit}",
                       initial=len(completed_tasks),
                       unit="tasks") as pbar:
            with ThreadPoolExecutor() as executor:
                future_to_task = {
                    executor.submit(
                        lambda t: t.wait(download_std_on_completion=False,
                                         silent_mode=True), task):
                        task for task in running_tasks
                }

                for future in as_completed(future_to_task):
                    with logging_redirect_tqdm():
                        task = future_to_task[future]
                        status = future.result()

                        if status != TaskStatusCode.SUCCESS:
                            logging.info(
                                "Task %s completed with status: %s\n"
                                "   · To understand why the task did not "
                                "complete successfully go to "
                                "https://console.inductiva.ai/tasks/%s",
                                task.id, status, task.id)

                        # Update progress bar for each completed task
                        pbar.update(1)

        logging.info("\n")
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

        if isinstance(select, str):
            select = SelectMode[select.upper()]

        if status is not None:
            status = TaskStatusCode(status)

        response = self._api.get_tasks_info_without_preload_content(
            name=self.name,
            select=select,
            status=status,
        )
        info = json.loads(response.data)

        if not filename:
            filename = f"{self.name}.{fmt.value}"
        elif not filename.endswith(f".{fmt.value}"):
            filename = f"{filename}.{fmt.value}"

        if fmt == ExportFormat.JSON:
            with open(filename, mode="w", encoding="utf-8") as file:
                json_content = json.dumps(obj=info, indent=4)
                file.write(json_content)
        elif fmt == ExportFormat.CSV:
            with open(filename, mode="w", encoding="utf-8") as file:
                fieldnames = sorted(info[0].keys()) if info else []
                writer = csv.DictWriter(file, fieldnames=fieldnames)
                writer.writeheader()
                writer.writerows(info)
        else:
            raise ValueError(f"Unsupported export format: {fmt}")

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
