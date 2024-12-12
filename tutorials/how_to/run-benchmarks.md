# Run a Benchmark

Benchmarking helps you measure and compare the performance of different configurations, 
like machine types or simulation settings. It’s an essential step to reduce 
costs and speed up your simulations without losing accuracy.

In this tutorial, we’ll show you how to use Inductiva’s API to run a benchmark. 
We’ll use the [SPlisHSPlasH simulator](https://tutorials.inductiva.ai/simulators/SPlisHSPlasH.html) 
as an example, so you can follow along and learn the process step by step.

We’ll benchmark the SPlisHSPlasH simulator across different machine types 
available on Google Cloud Platform (GCP). Specifically, we’ll test the 
performance of the `c2-standard` and `c3-standard` machine families with varying 
numbers of vCPUs. We’ll also compare these results to other configurations, 
such as `n2-standard-32`.

The goal is to find the best machine that balances computation time and cost. 
After running the benchmarks, we’ll analyze and visualize the results to make an 
informed choice.

## What You’ll Learn

1. How to set up a benchmark using Inductiva’s API.
2. How to configure and execute benchmarks with different machine types.
3. How to analyze the results to optimize your simulations.

By the end, you’ll know how to use benchmarks to optimize your computational resources. 

Let’s dive in!

## Step 1: Download the Files

Before we start benchmarking, we need to download the necessary input files 
for the SPlisHSPlasH simulator. These files include the configuration and 
assets required to run the simulations.

The code snippet below downloads the required simulation files to your working 
directory: 

```python3
import inductiva

inductiva.utils.download_from_url(
    url="https://tutorials.inductiva.ai/_static/generating-synthetic-data/splishsplash-base-dir.zip",
    unzip=True)
```

Now you’re set to proceed with setting up the benchmark!

## Step 2: Configure and run the benchmark (```run.py```)

- When initializing a benchmark, assign a name to the benchmark (e.g., ```splishsplash-fluid-cube```).
- Add new runs by calling ```add_run```.
- Execute the benchmark by calling ```run```, repeating each added run twice (i.e., ```num_repeats=2```).

```python
from inductiva import benchmarks, simulators, resources

benchmarks.Benchmark(name="splishsplash-fluid-cube") \
    .add_run(simulator=simulators.SplishSplash(),
             input_dir="splishsplash-base-dir",
             sim_config_filename="config.json",
             on=resources.MachineGroup("c2-standard-4")) \
    .add_run(simulator=simulators.SplishSplash(),
             input_dir="splishsplash-base-dir",             
             sim_config_filename="config.json",
             on=resources.MachineGroup("c2-standard-8")) \
    .add_run(simulator=simulators.SplishSplash(),
             input_dir="splishsplash-base-dir",             
             sim_config_filename="config.json",
             on=resources.MachineGroup("c2-standard-16")) \
    .add_run(simulator=simulators.SplishSplash(),
             input_dir="splishsplash-base-dir",
             sim_config_filename="config.json",
             on=resources.MachineGroup("c2-standard-30")) \
    .add_run(simulator=simulators.SplishSplash(),
             input_dir="splishsplash-base-dir",
             sim_config_filename="config.json",
             on=resources.MachineGroup("c2-standard-60")) \
    .add_run(simulator=simulators.SplishSplash(),
             input_dir="splishsplash-base-dir",
             sim_config_filename="config.json",
             on=resources.MachineGroup("c3-standard-4")) \
    .add_run(simulator=simulators.SplishSplash(),
             input_dir="splishsplash-base-dir",
             sim_config_filename="config.json",
             on=resources.MachineGroup("c3-standard-8")) \
    .add_run(simulator=simulators.SplishSplash(),
             input_dir="splishsplash-base-dir",
             sim_config_filename="config.json",
             on=resources.MachineGroup("c3-standard-22")) \
    .add_run(simulator=simulators.SplishSplash(),
             input_dir="splishsplash-base-dir",
             sim_config_filename="config.json",
             on=resources.MachineGroup("c3-standard-44")) \
    .add_run(simulator=simulators.SplishSplash(),
             input_dir="splishsplash-base-dir",
             sim_config_filename="config.json",
             on=resources.MachineGroup("c3-standard-88")) \
    .add_run(simulator=simulators.SplishSplash(),
             input_dir="splishsplash-base-dir",
             sim_config_filename="config.json",
             on=resources.MachineGroup("c3-standard-176")) \
    .run(num_repeats=2)
```

### Enhance the **readability** of the benchmark program

- Reduce lines of code, improve readability, and avoid code duplication by calling ```set_default```.

```python
from inductiva import benchmarks, simulators, resources

benchmarks.Benchmark(name="splishsplash-fluid-cube") \
    .set_default(simulator=simulators.SplishSplash(),
                 input_dir="splishsplash-base-dir",
                 sim_config_filename="config.json") \
    .add_run(on=resources.MachineGroup("c2-standard-4")) \
    .add_run(on=resources.MachineGroup("c2-standard-8")) \
    .add_run(on=resources.MachineGroup("c2-standard-16")) \
    .add_run(on=resources.MachineGroup("c2-standard-30")) \
    .add_run(on=resources.MachineGroup("c2-standard-60")) \
    .add_run(on=resources.MachineGroup("c3-standard-4")) \
    .add_run(on=resources.MachineGroup("c3-standard-8")) \
    .add_run(on=resources.MachineGroup("c3-standard-22")) \
    .add_run(on=resources.MachineGroup("c3-standard-44")) \
    .add_run(on=resources.MachineGroup("c3-standard-88")) \
    .add_run(on=resources.MachineGroup("c3-standard-176")) \
    .run(num_repeats=2)
```

- Create a list of machine types, then use a ```for``` loop to initialize the machine groups and add the new runs to the benchmark.

```python
from inductiva import benchmarks, simulators, resources

benchmark = benchmarks.Benchmark(name="splishsplash-fluid-cube") \
    .set_default(simulator=simulators.SplishSplash(),
                 input_dir="splishsplash-base-dir",
                 sim_config_filename="config.json")

machine_types = ["c2-standard-4", "c2-standard-8", "c2-standard-16",
                 "c2-standard-30", "c2-standard-60", "c3-standard-4",
                 "c3-standard-8", "c3-standard-22", "c3-standard-44",
                 "c3-standard-88", "c3-standard-176"]

for machine_type in machine_types:
    benchmark.add_run(on=resources.MachineGroup(machine_type))

benchmark.run(num_repeats=2)
```

### Optimize the benchmark program to minimize computation time and cost

- Avoid uploading the input files to each run by uploading them only once using the ```remote_assets``` parameter. This way, the input files are stored in a GCP bucket and reused on every simulation added to the benchmark.
- Begin by uploading the input files to a GCP bucket.

```python
import inductiva

inductiva.storage.upload(local_path="splishsplash-base-dir",
                         remote_dir="splishsplash-input-dir")
```

- Next, reuse the uploaded files on each run by passing the ```remote_assets``` argument to the ```set_default``` method.
- Remove the ```input_dir``` parameter from the ```set_default``` method.

```python
from inductiva import benchmarks, simulators, resources

benchmark = benchmarks.Benchmark(name="splishsplash-fluid-cube") \
    .set_default(simulator=simulators.SplishSplash(),
                 sim_config_filename="config.json",
                 remote_assets=["splishsplash-input-dir"])

machine_types = ["c2-standard-4", "c2-standard-8", "c2-standard-16",
                 "c2-standard-30", "c2-standard-60", "c3-standard-4",
                 "c3-standard-8", "c3-standard-22", "c3-standard-44",
                 "c3-standard-88", "c3-standard-176"]

for machine_type in machine_types:
    benchmark.add_run(on=resources.MachineGroup(machine_type))

benchmark.run(num_repeats=2)
```

- Speed up the benchmark execution by parallelizing the repetitions on each machine.
- For each machine group, set the number of machines to the number of repetitions (i.e., ```num_machines=num_repeats```)
- The repetitions of the simulations on each machine group will then run in parallel.

```python
from inductiva import benchmarks, simulators, resources

benchmark = benchmarks.Benchmark(name="splishsplash-fluid-cube") \
    .set_default(simulator=simulators.SplishSplash(),
                 sim_config_filename="config.json",
                 remote_assets=["splishsplash-input-dir"])

machine_types = ["c2-standard-4", "c2-standard-8", "c2-standard-16",
                 "c2-standard-30", "c2-standard-60", "c3-standard-4",
                 "c3-standard-8", "c3-standard-22", "c3-standard-44",
                 "c3-standard-88", "c3-standard-176"]

num_repeats = 2

for machine_type in machine_types:
    benchmark.add_run(on=resources.MachineGroup(machine_type=machine_type,
                                                num_machines=num_repeats))

benchmark.run(num_repeats=num_repeats)
```

- Reduce the maximum idle time of each resource to avoid wasting computational resources and, as a result, decrease the benchmark cost (i.e., ```max_idle_time = datetime.timedelta(seconds=30)```).
- Be aware that decreasing the maximum idle time may cause an error if there isn't enough time to submit the simulations.

```python
import datetime
from inductiva import benchmarks, simulators, resources

benchmark = benchmarks.Benchmark(name="splishsplash-fluid-cube") \
    .set_default(simulator=simulators.SplishSplash(),
                 sim_config_filename="config.json",
                 remote_assets=["splishsplash-input-dir"])

machine_types = ["c2-standard-4", "c2-standard-8", "c2-standard-16",
                 "c2-standard-30", "c2-standard-60", "c3-standard-4",
                 "c3-standard-8", "c3-standard-22", "c3-standard-44",
                 "c3-standard-88", "c3-standard-176"]

num_repeats = 2

max_idle_time = datetime.timedelta(seconds=30)

for machine_type in machine_types:
    benchmark.add_run(on=resources.MachineGroup(machine_type=machine_type,
                                                num_machines=num_repeats,
                                                max_idle_time=max_idle_time))

benchmark.run(num_repeats=num_repeats)
```

## Step 3: Export the benchmark data to a file (```export.py```)

- Export the metrics for the benchmark called ```"splishsplash-fluid-cube"``` into a CSV file named ```"benchmark.csv"```, filtering for only the successful tasks (```status="success"```) and selecting only the attributes and metrics that vary with each run (```select="distinct"```).
- To ensure that all benchmark-related tasks are completed and all resources used during the benchmark are terminated, call ```wait()``` and ```terminate()```, in that order.

```python
from inductiva import benchmarks

benchmarks.Benchmark(name="splishsplash-fluid-cube") \
    .wait() \
    .terminate() \
    .export(fmt="csv",
            status="success",
            select="distinct",
            filename="benchmark.csv")
```

## Step 4: Choose a visualization tool to analyze the benchmark results

- Open [CSVPlot](https://www.csvplot.com/) and upload the exported CSV file containing the benchmark data.

![plot](../_static/how_to/plot-benchmark-results.png)

- The best machine in terms of both computation time and cost is the ```c3-standard-44```, located at the bottom-left of the plot. This machine takes around 35 seconds and costs about $0.01, while the second-fastest machine, ```c3-standard-88```, takes the same amount of time but costs twice as much.

## Extra: Add another GCP instance to the benchmark for comparison with the previously added instances

- The user can add new runs later if desired.
- Recreate the benchmark using the same name (i.e., ```splishsplash-fluid-cube```).
- Add a new run execution for the ```n2-standard-32``` machine type (i.e., ```machine_type="n2-standard-32"```) and execute only this run twice.

```python
import datetime
from inductiva import benchmarks, simulators, resources

num_repeats = 2

max_idle_time = datetime.timedelta(seconds=30)

benchmark = benchmarks.Benchmark(name="splishsplash-fluid-cube") \
    .add_run(simulator=simulators.SplishSplash(),
             sim_config_filename="config.json",
             remote_assets=["splishsplash-input-dir"],
             on=resources.MachineGroup(machine_type="n2-standard-32",
                                       num_machines=num_repeats,
                                       max_idle_time=max_idle_time)) \
    .run(num_repeats=num_repeats)
```