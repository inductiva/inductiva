# Benchmarks

Our **Benchmark** tool is designed to help you run and evaluate simulations systematically, 
measuring their performance, cost, and execution time. It allows you to configure 
benchmarking runs, export results, and manage resources. 

## Quick Start
The typical workflow for running benchmarks follows three simple steps:
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

## Understanding Benchmarks

### Why Benchmark?
Benchmarking is essential for making informed decisions about computational resources. It helps you identify the optimal balance between performance, cost, and accuracy for your specific simulation requirements. Without systematic benchmarking, you might overpay for unnecessary compute power or accept poor performance from inadequate resources.

### When to Benchmark
Consider benchmarking when you need to:
- Choose between different machine types
- Evaluate the impact of simulation parameters on performance
- Understand cost-performance trade-offs across different configurations
- Validate that performance improvements justify increased costs

### How Benchmarks Work
The `Benchmark` API class handles the complexity of resource provisioning, parallel execution, and result aggregation, allowing you to focus on defining meaningful test scenarios.

## `Benchmark` Class
### Core Concepts

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

````{eval-rst}
.. seealso::
   For complete API documentation including all parameters, methods, and configuration options, see the `Benchmark <https://inductiva.ai/guides/api-functions/api/inductiva.benchmarks>`_ class documentation
````

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

### Design Principles

**Start Simple, Then Expand**
Begin with a small set of representative configurations before scaling to comprehensive evaluations. This helps validate your approach and avoid costly mistakes.

**Control Variables**
Change only one parameter at a time when possible. This makes it easier to identify the source of performance differences and builds clearer insights.

**Representative Workloads**
Use simulation inputs that represent your actual production workloads. Synthetic benchmarks might not capture real-world performance characteristics.

### Resource Management

**Always Clean Up**
```python
# Ensure resources are terminated after benchmarking
benchmark.wait().terminate()
```

**Monitor Costs**

Regularly review resource usage to avoid unexpected charges during long-running benchmarks.

**Use Quotas Wisely**

Consider enabling `wait_for_quotas=True` in `run()` to ensure your simulation runs without quota issues.
```python
# Handle resource quotas gracefully
benchmark.run(num_repeats=3, wait_for_quotas=True)
```

**Optimize Idle Time**
```python
# Minimize idle time to reduce costs
max_idle_time = datetime.timedelta(seconds=30)
```

### Statistical Considerations

**Account for Variability**

Run multiple repetitions to ensure reliable results. System variability can significantly impact single-run measurements:

```python
# Minimum recommended repetitions
benchmark.run(num_repeats=3)

# For critical benchmarks, use more repetitions
benchmark.run(num_repeats=5)
```

### Performance Optimization

**Parallel Repetitions**

Optimize benchmark execution time by running repetitions in parallel:

```python
num_repeats = 3

for machine_type in machine_types:
    benchmark.add_run(
        on=resources.MachineGroup(
            machine_type=machine_type,
            num_machines=num_repeats  # Each repetition runs on separate machine
        )
    )

benchmark.run(num_repeats=num_repeats)
```

---

Ready to dive in? Check out these exciting tutorials and our blog for more insights:

- [Allocating Computational Resources in a Diverse Chip Ecosystem](https://inductiva.ai/blog/article/allocating-computational-resources-in-a-diverse-chip-ecosystem)
- [How to Run a Benchmark](https://tutorials.inductiva.ai/how_to/run-benchmarks.html)
- [Benchmarking Computational Resources - A Use Case](https://tutorials.inductiva.ai/generating-synthetic-data/synthetic-data-generation-6.html)
