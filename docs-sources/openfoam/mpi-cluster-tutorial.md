# Scaling the MB9 Microbenchmark from ExaFOAM Using MPI Clusters
In this guide, you'll learn how to scale an advanced OpenFOAM case across **multiple machines** using an **MPI (Message Passing Interface) cluster**. This approach combines the CPU resources of several machines, significantly boosting the total number of virtual CPUs (vCPUs) available and accelerating your simulation’s runtime.

With **Inductiva**, running your simulations on a multi-node MPI cluster is as straightforward as 
using any single-node option.

All you need to do is update the resource allocation from `MachineGroup` to `MPICluster`, as follows:

```diff
# Allocate a multi-machine MPI cluster on Google Cloud Platform
- cloud_machine = inductiva.resources.MachineGroup(
+ cloud_machine = inductiva.resources.MPICluster(
    machine_type="c3d-highcpu-360",
    num_machines=3,
    data_disk_gb=300,
    spot=True
```

In this tutorial, we’ll demonstrate how leveraging an MPI cluster can impact the performance of your OpenFOAM simulation. To do this, we will revisit the [MB9 Microbenchmark from ExaFOAM](https://exafoam.eu/benchmarks/), which was previously explored in [this tutorial](run-exafoam-microbenchmark).

## Prerequisites

### Adjust the Case
Begin by completing the setup steps outlined in the Prerequesites chapter of our [tutorial](run-exafoam-microbenchmark).

Instead of 360, update the `nCores` value in `highLiftConfiguration/system/include/caseDefinition` to 1080, since it's the total number of vCPUs the MPI Cluster will have.

### Convert the `Allrun` Script into Python Commands
To fully leverage MPI clusters, you **cannot** use the `shell_script` argument as you might with OpenFOAM simulations. 
Instead, you must use the `commands` argument, which accepts a list of commands to be executed during the simulation. 
Some commands will run sequentially, while others may run in parallel.

For this case, we’ve prepared a complete list of required commands. Simply download the `commands.txt` [here](https://storage.googleapis.com/inductiva-api-demo-files/commands.txt) and place it in your `highLiftConfiguration` directory. 

If you’d prefer to do it yourself, check out our [documentation](https://inductiva.ai/guides/openfoam/convert-allrun-script-into-python-commands)

## Running the Simulation on an MPI Cluster
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

Estimated computation cost (US$): 179.23 US$
```

As you can see in the "In Progress" line, the part of the timeline that represents the actual execution of the simulation, 
the core computation time of this simulation was around **17 hours minutes and 11 minutes**.

## Performance Comparison
Below is a comparison of the same simulation run across different MPI cluster configurations:

| Machine Type    | Machines | Total Cores | Duration                 |
| --------------- | -------- | ----------- | ------------------------ |
| c3d-highcpu-360 | 1        | 360         | 41 h                     |
| c3d-highcpu-360 | 3        | 1080        | 17 h 11 min              |
| c2d-highcpu-112 | 8        | 896         | 14 h 59 min              |
| c2d-highcpu-112 | 12       | 1344        | 12 h 5 min *(Preempted)* |

As you can see, scaling up improves performance significantly. However, the gains diminish at higher core counts — going from 896 to 1344 cores didn't offer much improvement.



