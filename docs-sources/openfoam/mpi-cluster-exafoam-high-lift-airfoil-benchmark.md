# Scaling the 30P30N Micro-Benchmark (MB9) from ExaFOAM Using MPI Clusters

In this guide, you'll learn how to scale the
[30P30N micro-benchmark](run-exafoam-high-lift-airfoil-benchmark) from ExaFOAM
across multiple machines using MPI clusters.

MPI (Message Passing Interface) clusters allow you to link several machines into
one very powerful machine. For example, if you provision an MPI cluster with 2
machines, each with 100 cores, you can run a simulation using all 200 cores as
if they belonged to a single machine.

## 🚧 Prerequisites

Start by completing the setup steps from the [original tutorial](run-exafoam-high-lift-airfoil-benchmark).

To take full advantage of MPI clusters, you **cannot** use the `shell_script`
argument as before. Instead, you'll use the `commands` argument, which takes a
list of commands that will be executed during the simulation. Some will
run sequentially, others in parallel.

## 🔁 Converting the `Allrun` Script into Python Commands

Transforming your `Allrun` script into a list of commands is simple. Just
extract each command from the script (without logic or flow control) and insert
it into a Python list as strings.

### Example Conversion

The following shell snippet:

```bash
echo -e -n '\n    Running snappyHexMesh ...\n'
(
    cd system
    cp controlDict.SHM controlDict
    cp fvSchemes.SHM fvSchemes
    cp fvSolution.SHM fvSolution
    cd ..
)
```

Becomes this in Python:

```python
commands = [
    # You can also add the echo if you wanted
    "cp system/controlDict.SHM system/controlDict",
    "cp system/fvSchemes.SHM system/fvSchemes",
    "cp system/fvSolution.SHM system/fvSolution",
]
```

> ⚠️ **Important Notes**
>
> * Each command must be independent. You **cannot** use `cd` to change directories between commands.
> * Logic such as conditionals or loops must be removed.
> * Special characters like `>`, `|`, and `&` are **not allowed**.

We've prepared a complete list of required commands for you.
👉 [Download commands.txt](https://storage.googleapis.com/inductiva-api-demo-files/commands.txt) and place it in your `highLiftConfiguration` directory.

Finally, update the `nCores` value in `highLiftConfiguration/system/include/caseDefinition` to 1080 (1080 is the total number of vcpus your MPI Cluster will have).

## 🚀 Running the Simulation on an MPI Cluster

Here's how to launch your simulation using the Inductiva API:

```python
"""Run OpenFOAM on an MPI Cluster"""
import inductiva

# Set up an MPI Cluster on Google Cloud
cloud_machine = inductiva.resources.MPICluster(
    machine_type="c3d-highcpu-360",
    num_machines=3,
    data_disk_gb=300,
    spot=True
)

# Initialize the simulator
OpenFOAM = inductiva.simulators.OpenFOAM(
    version="2412",
    distribution="esi"
)

# Load command list
with open("/path/to/highLiftConfiguration/commands.txt", "r") as file:
    commands = [line.strip() for line in file]

# Run the simulation
task = OpenFOAM.run(
    input_dir="/path/to/highLiftConfiguration",
    commands=commands,
    on=cloud_machine
)

# Wait for completion and fetch results
task.wait()
cloud_machine.terminate()
task.download_outputs()
task.print_summary()
```



## 📊 Task Summary

When the simulation completes, you'll get a summary like this:

```
inductiva tasks info hvtitxq0h6aw9yagbfhdwsc5v

Task status: Success

Timeline:
	Waiting for Input         at 23/04, 16:46:16      2.615 s
	In Queue                  at 23/04, 16:46:19      89.779 s
	Preparing to Compute      at 23/04, 16:47:48      2.87 s
	In Progress               at 23/04, 16:47:51      61841.482 s
        ...
    Finalizing                at 24/04, 09:58:33      817.51 s
	Success                   at 24/04, 10:12:10      

Data:
	Size of zipped output:    26.72 GB
	Size of unzipped output:  34.41 GB
	Number of output files:   151512

Estimated computation cost (US$): 179.23 US$

Go to https://console.inductiva.ai/tasks/hvtitxq0h6aw9yagbfhdwsc5v for more details.
```

The "In Progress" section shows the actual simulation time — in this example, around **17 hours and 11 minutes**.

## ⚖️ Performance Comparison

Below is a comparison of the same simulation run across different MPI cluster configurations:

| Machine Type    | Machines | Total Cores | Duration                 |
| --------------- | -------- | ----------- | ------------------------ |
| c3d-highcpu-360 | 1        | 360         | 41 h                     |
| c3d-highcpu-360 | 3        | 1080        | 17 h 11 min              |
| c2d-highcpu-112 | 8        | 896         | 14 h 59 min              |
| c2d-highcpu-112 | 12       | 1344        | 12 h 5 min *(Preempted)* |

As you can see, scaling up improves performance significantly. However, the gains diminish at higher core counts — going from 896 to 1344 cores didn't offer much improvement.

## 🧪 Final Notes

Feel free to experiment with different cluster configurations to find
the sweet spot between performance and cost. Good luck, and happy simulating! 🎉


