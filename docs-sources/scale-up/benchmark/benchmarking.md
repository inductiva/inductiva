# Benchmarks

The Inductiva API provides a **benchmarking tool** to help you measure, compare, and optimize simulation performance across different configurations. Make data-driven decisions about resource allocation, cost optimization, and performance tuning.

## Quick Start
The `Benchmark` class is the core of this tool, handling the complexity of resource provisioning, managing runs, and aggregating results, so you can focus on defining your test scenarios. A typical workflow involves three simple steps:
1. **Setup** your benchmark configuration with different runs to compare
2. **Execute** the benchmark runs and monitor their progress
3. **Analyze** results and terminate resources to avoid unnecessary costs

```python
from inductiva import simulators, resources, benchmarks

# Step 1: Setup - Create and configure your benchmark
benchmark = benchmarks.Benchmark(name="MyFirstBenchmark")
benchmark.set_default(
    simulator=simulators.OpenFOAM(),
    input_dir="/path/to/input",
    sim_config_filename="config_file"
)
benchmark.add_run(on=resources.MachineGroup("c2-standard-4"))
benchmark.add_run(on=resources.MachineGroup("c2-standard-8"))

# Step 2: Execute - Run the benchmark and monitor progress
benchmark.run(num_repeats=3)
benchmark.wait()

# Step 3: Analyze - Export results and clean up resources
benchmark.export(fmt="csv", filename="results.csv")
benchmark.terminate()
```

## `Benchmark` Class

A `Benchmark` object contains **runs**, where each run represents a different simulation configuration you want to compare. This might include different:
- Machine types or sizes
- Geographical zones
- Simulation parameters

### Key Methods

| Method | Purpose | When to Use |
|--------|---------|-------------|
| `set_default()` | Configure common parameters for all runs | Set simulator, base machine type, common inputs |
| `add_run()` | Add a specific configuration to test | Each variation you want to benchmark |
| `run()` | Execute all configured runs | After setting up all test scenarios |
| `wait()` | Block until all tasks complete | Before analyzing results |
| `export()` | Save results to file | For analysis and reporting |
| `terminate()` | Clean up resources | Always call when finished |


> For complete API documentation including all parameters, methods, and configuration options, see the [Benchmark](https://inductiva.ai/guides/api-functions/api/inductiva.benchmarks) class documentation

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

> 💡 Tip: Use `set_default()` to configure common parameters (e.g., simulator, resources) for all runs.

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
Always call `terminate()` after completing a benchmark to free up resources, especially if using cloud providers like GCP.

```py
# Wait for tasks to finish
benchmark.wait()

# Terminate all machine groups used during the benchmark
benchmark.terminate()
```

### Example 4: Exporting Benchmark Results

```py
# Export the benchmark results in CSV format
benchmark.export(fmt="csv", filename="benchmark_results.csv")

# Export only distinct data for successfull tasks in JSON format
benchmark.export(fmt="json", select="distinct", status="success")
```

> 💡 Tip: Export results regularly! After completing benchmark runs, export the results in your preferred format to facilitate analysis and reporting. 

## Best Practices

**Start Simple, Then Expand**
Begin with a small set of representative configurations before scaling to comprehensive evaluations. This helps validate your approach and avoid costly mistakes.

**Control Variables**
Change only one parameter at a time when possible. This makes it easier to identify the source of performance differences and builds clearer insights.

**Representative Workloads**
Use simulation inputs that represent your actual production workloads. Synthetic benchmarks might not capture real-world performance characteristics.

## Next Steps
Now that you understand the fundamentals of benchmarks in the Inductiva API, explore these topics to deepen your knowledge:

- [How to Run a Benchmark](run-benchmarks.md)
- [Allocating Computational Resources in a Diverse Chip Ecosystem](https://inductiva.ai/blog/article/allocating-computational-resources-in-a-diverse-chip-ecosystem)