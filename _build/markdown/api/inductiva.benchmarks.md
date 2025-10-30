# benchmarks

API for Benchmarking

### *class* Benchmark(name: str, verbose: bool = False)

Bases: [`Project`](inductiva.projects.md#inductiva.projects.project.Project)

Represents the benchmark runner.

#### \_\_init_\_(name: str, verbose: bool = False)

Initializes a new Benchmark instance.

* **Parameters:**
  **name** (*str*) – The name of the benchmark runner.

#### add_run(simulator: Simulator | None = None, input_dir: str | None = None, on: [MachineGroup](inductiva.resources.md#inductiva.resources.machine_groups.MachineGroup) | [ElasticMachineGroup](inductiva.resources.md#inductiva.resources.machine_groups.ElasticMachineGroup) | [MPICluster](inductiva.resources.md#inductiva.resources.machine_groups.MPICluster) | None = None, \*\*kwargs) → Self

Adds a simulation run to the benchmark.

* **Parameters:**
  * **simulator** (*Optional* *[**Simulator* *]*) – The simulator to be used for this
    run. If not provided, the previously set simulator will be used.
  * **input_dir** (*Optional* *[**str* *]*) – The directory containing input files for
    the simulation. If not provided, the previously set input
    directory will be used.
  * **on** (*Optional* *[**types.ComputationalResources* *]*) – The computational
    resources to run the simulation on. If not provided, the
    previously set resources will be used.
  * **\*\*kwargs** – Additional keyword arguments for the simulator run. These
    will overwrite any previously set parameters.
* **Returns:**
  The current instance for method chaining.
* **Return type:**
  Self

#### add_task(task: [Task](inductiva.tasks.md#inductiva.tasks.task.Task))

Adds a task to the project.

* **Parameters:**
  **task** – The task to add to the project.

#### *property* created_at *: str*

Returns the creation date and time of the project.

#### delete()

Delete a project on the backend.

This method does not delete the project tasks, only the project itself.
The tasks will be moved to the “default” project.

#### download_outputs(output_dir: str | None = None)

Downloads all the outputs for all the tasks in the project.

All task outputs will be organized within the specified output_dir.
If output_dir is not provided, outputs will be saved to a default
location under inductiva_output/<project_name>/<task_id>/.
Otherwise, they will be stored in <output_dir>/<task_id>/.

* **Parameters:**
  **output_dir** (*str* *,* *optional*) – The base directory where project outputs
  will be downloaded.

#### *property* estimated_computation_cost *: float*

Returns the estimated project cost.

Computed as the sum of the estimated computation cost of each task.

#### export(fmt: [ExportFormat](#inductiva.benchmarks.benchmark.ExportFormat) | str = ExportFormat.JSON, filename: str | None = None, status: TaskStatusCode | str | None = None, select: [SelectMode](#inductiva.benchmarks.benchmark.SelectMode) | str = SelectMode.DISTINCT)

Exports the benchmark performance metrics in the specified format.

* **Parameters:**
  * **fmt** (*Union* *[*[*ExportFormat*](#inductiva.benchmarks.benchmark.ExportFormat) *,* *str* *]*) – The format to export the results
    in. Defaults to ExportFormat.JSON.
  * **filename** (*Optional* *[**str* *]*) – The name of the output file to save the
    exported results. Defaults to the benchmark’s name if not
    provided.
  * **status** (*Optional* *[**Union* *[**TaskStatusCode* *,* *str* *]* *]*) – The status of the
    tasks to include in the benchmarking results. Defaults to None,
    which includes all tasks.
  * **select** (*Union* *[*[*SelectMode*](#inductiva.benchmarks.benchmark.SelectMode) *,* *str* *]*) – The data to include in
    the benchmarking results. Defaults to SelectMode.DISTINCT that
    includes only the parameters that vary between different runs.

#### get_tasks(last_n: int = -1, status: str | None = None) → List[[Task](inductiva.tasks.md#inductiva.tasks.task.Task)]

Get the the tasks of this project.

Optionally, those can be filtered by task status.

* **Parameters:**
  * **last_n** (*int*) – The number of tasks with repect to the submission
    time to fetch. If last_n<=0 we fetch all tasks submitted
    to the project.
  * **status** – Status of the tasks to get. If None, tasks with any
    status will be returned.

#### *property* id *: str*

Returns the unique ID of the project.

#### *property* name *: str*

Returns the name of the project.

#### *property* num_tasks *: int*

Returns the number of tasks in the project.

#### run(num_repeats: int = 2, wait_for_quotas: bool = True) → Self

Executes all added runs.

Each run is executed the specified number of times, and the collection
of runs is cleared afterwards.

* **Parameters:**
  * **num_repeats** (*int*) – The number of times to repeat each simulation
    run (default is 2).
  * **wait_for_quotas** (*bool*) – Indicates whether to wait for quotas to
    become available before starting each resource. If True, the
    program will actively wait in a loop, periodically sleeping and
    checking for quotas. If False, the program crashes if quotas
    are not available (default is True).
* **Returns:**
  The current instance for method chaining.
* **Return type:**
  Self

#### set_default(simulator: Simulator | None = None, input_dir: str | None = None, on: [MachineGroup](inductiva.resources.md#inductiva.resources.machine_groups.MachineGroup) | [ElasticMachineGroup](inductiva.resources.md#inductiva.resources.machine_groups.ElasticMachineGroup) | [MPICluster](inductiva.resources.md#inductiva.resources.machine_groups.MPICluster) | None = None, \*\*kwargs) → Self

Sets default parameters for the benchmark runner.

This method allows you to configure default settings for the benchmark,
which will be used in subsequent runs unless explicitly overridden.

* **Parameters:**
  * **simulator** (*Optional* *[**Simulator* *]*) – The simulator instance to be used
    as the default for future runs. If not provided, the current
    simulator will remain unchanged.
  * **input_dir** (*Optional* *[**str* *]*) – The directory path for input files. If
    not provided, the current input directory will remain unchanged.
  * **on** (*Optional* *[**types.ComputationalResources* *]*) – The computational
    resources to use for running the simulations. If not specified,
    the current resources will remain unchanged.
  * **\*\*kwargs** – Additional keyword arguments to set as default parameters
    for the simulations. These will update any existing parameters
    with the same names.
* **Returns:**
  The current instance for method chaining.
* **Return type:**
  Self

#### *property* task_by_status *: dict*

Returns a dictionary with the number of tasks by status.
The keys are the status codes and the values are the number of tasks
with that status.

#### terminate() → Self

Terminates all active machine groups associated with the
benchmark.

* **Returns:**
  The current instance for method chaining.
* **Return type:**
  Self

#### *property* total_estimated_cost *: float*

#### *property* total_task_orchestration_fee *: float*

#### wait() → Self

Waits for all running tasks to complete.

* **Returns:**
  The current instance for method chaining.
* **Return type:**
  Self

### *class* ExportFormat(value)

Bases: `Enum`

Enumeration of supported benchmark export formats.

#### CSV *= 'csv'*

#### JSON *= 'json'*

### *class* SelectMode(value)

Bases: `Enum`

Enumeration of supported data selection modes, specifying which data
should be included in the benchmarking results.

#### ALL *= 'all'*

#### DISTINCT *= 'distinct'*

::docsbannersmall
::
