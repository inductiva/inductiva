# Benchmarking

The **Benchmarking API** is designed to help users run and evaluate simulations, 
measuring their performance, cost, and execution time. It allows you to configure 
benchmarking runs, export results, and manage resources. 

The tool supports different export formats, handles the execution of multiple 
simulations in parallel, and provides mechanisms to gather and present performance 
metrics.

Here you can find documentation on the benchmarking features of our Inductiva API 
in Python.

## Key Classes and Concepts

### **`Benchmark`**

The `Benchmark` class is the core of the benchmarking tool. It manages benchmarking 
runs, including specifying simulation parameters, running simulations, and exporting results.

#### Constructor

```py
Benchmark(name: str, append: bool = True)
```

* **name (str)**: The name of the benchmark. This will be used for identification and in output filenames.  
* **append (bool)**: Indicates whether to allow adding runs to the existing benchmark (default is `True`).

#### Methods

##### `set_default`

```py
set_default(simulator: Optional[simulators.Simulator] = None, 
            input_dir: Optional[str] = None, 
            on: Optional[types.ComputationalResources] = None, 
            **kwargs) -> Self
```

Sets default parameters for all benchmarking runs.

* **simulator (Optional\[simulators.Simulator\])**: The simulator to be used for future runs.  
* **input\_dir (Optional\[str\])**: Directory for input files for the simulation in future runs.  
* **on (Optional\[types.ComputationalResources\])**: Computational resources for running the simulation in future runs.  
* **kwargs**: Additional parameters to be passed to each run (e.g., simulation-specific settings).

##### **`add_run`**

```py
add_run(simulator: Optional[simulators.Simulator] = None, 
        input_dir: Optional[str] = None, 
        on: Optional[types.ComputationalResources] = None, 
        **kwargs) -> Self
```

Adds a new run to the benchmark. This allows you to specify the simulation parameters for each individual run.

* **simulator (Optional\[simulators.Simulator\])**: The simulator for the current run.  
* **input\_dir (Optional\[str\])**: Directory for input files for the run.  
* **on (Optional\[types.ComputationalResources\])**: Computational resources for the run.  
* **kwargs**: Additional parameters specific to this run.

##### **`run`**

```py
run(num_repeats: int = 2, wait_for_quotas: bool = True) -> Self
```

Runs all the added benchmarking simulations. Each simulation can be repeated multiple times.

* **num\_repeats (int)**: The number of times each simulation run should be repeated (default is 2).  
* **wait\_for\_quotas (bool)**: Whether to wait for resource quotas before starting the runs (default is `False`).

##### **`wait`**

```py
wait() -> Self
```

Waits for all running tasks to finish. This is useful when you want to ensure that all tasks have completed before proceeding.

##### **`export`**

```py
export(fmt: Union[ExportFormat, str] = ExportFormat.JSON, 
       filename: Optional[str] = None, 
       status: Optional[Union[TaskStatusCode, str]] = None, 
       select: Union[SelectMode, str] = SelectMode.DISTINCT)
```

Exports the benchmark performance data to a specified file in the desired format.

* **fmt (Union\[ExportFormat, str\])**: The format for the exported data. Can be `JSON` (default) or `CSV`.  
* **filename (Optional\[str\])**: The filename to save the exported data. If not provided, the default name based on the benchmark's name will be used.  
* **status (Optional\[Union\[TaskStatusCode, str\]\])**: The status of the tasks to be included in the export. If not provided, all tasks will be included.  
* **select (Union\[SelectMode, str\])**: Specifies which data to include in the export. `DISTINCT` (default) includes only distinct parameters that vary between runs, and `ALL` includes all parameters.

##### **`terminate`**

```py
terminate() -> Self
```

Terminates all active machine groups associated with the benchmark.

---

## Available Enums

### **`ExportFormat`**

Specifies the format for exporting benchmark results.

```py
class ExportFormat(enum.Enum):
    JSON = "json"
    CSV = "csv"
```

* `ExportFormat.JSON`: Exports results in JSON format.  
* `ExportFormat.CSV`: Exports results in CSV format.

### **`SelectMode`**

Specifies the selection mode for the exported benchmark data.

```py
class SelectMode(enum.Enum):
    ALL = "all"
    DISTINCT = "distinct"
```

* `SelectMode.ALL`: Exports all data for the runs.  
* `SelectMode.DISTINCT`: Exports only distinct parameters that vary between different runs.

---

## Code Examples

### Example 1: Creating a Basic Benchmark

```py
from inductiva import simulators, resources, benchmarks

# Create a benchmark
benchmark = benchmarks.Benchmark(name="MyBenchmark")

# Set default simulation parameters
benchmark.set_default(
    simulator=simulators.OpenFOAM(),
    input_dir="/path/to/input",
    on=resources.MachineGroup("c2-standard-4"),
)

# Add a new run with specific parameters
benchmark.add_run(input_dir="/path/to/another/input",
                  num_iterations=100)

# Run the benchmark
benchmark.run(num_repeats=3)
```

### Example 2: Adding Multiple Runs with Different Parameters

```py
# Add multiple runs with different simulators and parameters
benchmark.add_run(simulator=simulators.AnotherSimulator(),
                  input_dir="/input/dir1",
                  on=resources.MachineGroup("c2-standard-8"))
benchmark.add_run(simulator=simulators.MySimulator(),
                  input_dir="/input/dir2",
                  on=resources.MachineGroup("c2-standard-16"),
                  num_iterations=200)

# Execute all runs twice
benchmark.run(num_repeats=2)
```

### Example 3: Waiting for Tasks and Terminating Resources

```py
# Wait for tasks to finish
benchmark.wait()

# Terminate all machine groups used during the benchmark
benchmark.terminate()
```

### Example 4: Exporting Benchmark Results

```py
# Export the benchmark results in CSV format
benchmark.export(
    fmt="csv",
    filename="benchmark_results.csv")

# Export only distinct data for successfull tasks in JSON format
benchmark.export(
    fmt="json",
    select="distinct",
    status="success")
```

---

**Best Practices**

* **Set Default Parameters**: Use `set_default()` to configure common parameters (e.g., simulator, resources) for all runs.  
* **Export Results Regularly**: After completing benchmark runs, export the results in your preferred format to facilitate analysis and reporting.  
* **Terminate Resources**: Always call `terminate()` after completing a benchmark to free up resources, especially if using cloud providers like GCP.  
* **Handle Quotas**: Consider enabling `wait_for_quotas=True` in `run()` to ensure your simulation runs without quota issues.

---

Ready to dive in? Check out these exciting tutorials and our blog for more insights:

- [Allocating Computational Resources in a Diverse Chip Ecosystem](https://inductiva.ai/blog/article/allocating-computational-resources-in-a-diverse-chip-ecosystem)
- [How to Run a Benchmark](https://tutorials.inductiva.ai/how_to/run-benchmarks.html)
- [Benchmarking Computational Resources - A Use Case](https://tutorials.inductiva.ai/generating-synthetic-data/synthetic-data-generation-6.html)
