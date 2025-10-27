# Run the MB9 Microbenchmark from ExaFOAM on an MPI Cluster
In this tutorial, you’ll learn how to run a high-performance OpenFOAM simulation on a **multi-node MPI (Message Passing Interface) cluster** using the Inductiva API. This setup lets you combine the CPU resources of multiple machines, significantly increasing the total number of virtual CPUs (vCPUs) and reducing simulation runtime.

As a practical example, we’ll run the `MB9 Microbenchmark from ExaFOAM`, available in the official [official ExaFOAM documentation](https://exafoam.eu/benchmarks/), using three machines within the cluster. This configuration **cuts the simulation time in half** compared to using a single machine.

## Prerequisites
1. Download the required files [here](https://develop.openfoam.com/committees/hpc/-/tree/develop/compressible/rhoPimpleFoam/LES/highLiftConfiguration) and place them in a folder named `highLiftConfiguration`.

The **directory structure** should look like this:
```bash
ls -lasgo highLiftConfiguration
total 128
 0 drwxrwxr-x@ 13     416 Apr 10 11:13 .
 0 drwx------@  5     160 Apr 23 08:35 ..
24 -rw-r--r--@  1   10244 Jun 12 15:58 .DS_Store
 0 drwxrwxr-x@ 11     352 Apr  8 11:27 0.orig
 8 -rwxr-xr-x@  1     626 Apr  8 11:27 Allclean
16 -rwxr-xr-x@  1    7017 Jun 12 15:41 Allrun
 8 -rw-rw-r--@  1     991 Apr  8 11:27 COPYING
48 -rw-rw-r--@  1   21547 Apr  8 11:27 README.md
 0 -rw-rw-r--@  1       0 Apr  8 11:27 case.foam
 0 drwxrwxr-x@  5     160 Apr  8 11:27 constant
 0 drwxrwxr-x@ 13     416 Apr  8 11:27 figures
 0 drwxrwxr-x@ 28     896 Apr 10 11:13 system
24 -rw-rw-r--@  1   11399 Apr  8 11:27 thumbnail.png
```

2. Apply the following modifications to the `Allrun` and `caseDefinition` files.

* Open the `Allrun` script and **replace**:

```bash
parEx="mpirun -np $nProcs"
```

with:

```bash
parEx="mpirun -use-hwthread-cpus -np $nProcs"
```

> The `-use-hwthread-cpus` flag enables all available virtual CPUs on the machine for optimal performance.

* Edit `highLiftConfiguration/system/include/caseDefinition` and **update** the following parameters:
- Set Time Step (`dt`) to 0.00002
- Set Start Time (`initTime`) to 0.10
- Set End Time (`finalTime`) to 0.30
- Set Number of Cores (`nCores`) to 1080 - total number of vCPUs in your MPI cluster

These modifications are made for testing purposes, in order to shorten the
overall simulation time while still allowing you to validate the setup and
performance.

3. Convert the `Allrun` Script into Python Commands.

To fully leverage MPI clusters, you **cannot** use the `shell_script` argument
as you might with OpenFOAM simulations. This is because if we use a shell script,
Inductiva cannot see or analyze its contents, it just executes it as-is. For MPI
runs on a cluster, we need to prepare commands (e.g., `mpirun`) with additional
configurations so different machines can communicate, such as specifying a
machinefile with all cluster nodes.

Instead, you must use the `commands` argument, which accepts a list of commands
to be executed during the simulation. This way, Inductiva can process each
command and apply the necessary cluster-specific options. Some commands will
run sequentially, while others may run in parallel.

For this case, we’ve prepared a list of the required commands. You can download the `commands.txt` [here](https://storage.googleapis.com/inductiva-api-demo-files/commands.txt) and place it in your `highLiftConfiguration` directory.

Prefer to create it yourself? Check out our guide on [converting an Allrun script into Python commands](convert-allrun-script-into-python-commands).

## Running the Simulation on an MPI Cluster
With Inductiva, running simulations on a multi-node MPI cluster is just as simple as using a single-node setup.
The only change required is to switch the resource allocation from `MachineGroup` to `MPICluster`, as shown below:

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

# Initialize the Simulator
OpenFOAM = inductiva.simulators.OpenFOAM(
    version="2412",
    distribution="esi"
)

# Load command list
with open("/path/to/highLiftConfiguration/commands.txt", "r") as file:
    commands = [line.strip() for line in file]

# Run simulation
task = OpenFOAM.run(
    input_dir="/path/to/highLiftConfiguration",
    commands=commands,
    on=cloud_machine
)

# Wait for the simulation to finish and download the results
task.wait()
cloud_machine.terminate()
task.download_outputs()
task.print_summary()
```

We run the simulation specifying the command list that handles the execution process.

When the simulation is complete, we terminate the machine, download the results and print a summary of the simulation as shown below.

```
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

Total estimated cost (US$): 179.24 US$
	Estimated computation cost (US$): 179.23 US$
	Task orchestration fee (US$): 0.010 US$

Note: A per-run orchestration fee (0.010 US$) applies to tasks run from 01 Dec 2025, in addition to the computation costs.
Learn more about costs at: https://inductiva.ai/guides/how-it-works/basics/how-much-does-it-cost
```

As you can see in the "In Progress" line, the part of the timeline that represents the actual execution of the simulation, the core computation time of this simulation was around **17 hours and 11 minutes**.

## Performance Comparison
The table below compares the performance of the same simulation running on a single machine and on a three-machine MPI cluster configurations.

| Machine Type    | Nº of Machines | Total Cores | Execution Time           |
| --------------- | -------------- | ----------- | ------------------------ |
| c3d-highcpu-360 | 1              | 360         | 41h                      |
| c3d-highcpu-360 | 3              | 1080        | 17h, 11 min              |

As shown, using an MPI cluster with three machines reduced the execution time
from 41 hours to just over 17 hours, demonstrating the significant performance
benefits of leveraging multiple nodes for high-performance computing tasks.

⚠️ Keep in mind that MPI clusters introduce communication overhead between
machines. This means scaling is not perfectly linear, adding more machines
doesn’t always lead to proportional speedups, and in some cases, performance
may even degrade if the communication cost becomes too high.



