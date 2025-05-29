# Running AMR-Wind Simulations Across Multiple Machines Using MPI
As computational fluid dynamics (CFD) problems increase in complexity and scale, so
does the demand for computing power. **AMR-Wind** is designed to scale efficiently across thousands of cores, but even the most powerful single machine can become a bottleneck, particularly when handling large domains or very fine mesh resolutions.

To overcome these limitations, you can run AMR-Wind simulations across
**multiple machines** using an **MPI (Message Passing Interface) cluster**. This approach combines the 
CPU resources of several machines, significantly boosting the total number of virtual CPUs (vCPUs) 
available and accelerating your simulation’s runtime.

With **Inductiva**, running your simulations on a multi-node MPI cluster is as straightforward as 
using any single-node option.

All you need to do is update the resource allocation from `MachineGroup` to `MPICluster`, like this:

```diff
# Allocate a multi-machine MPI cluster on Google Cloud Platform
- cloud_machine = inductiva.resources.MachineGroup(
+ cloud_machine = inductiva.resources.MPICluster(
    machine_type="c2d-highcpu-112",
+   num_machines=2,
    data_disk_gb=100,
    spot=True)
```

In this tutorial, we’ll demonstrate how leveraging an MPI cluster can impact the performance of your AMR-Wind simulation 
by using the neutral Atmospheric Boundary Layer case from the [AMR-Wind GitHub repository](https://github.com/Exawind/amr-wind/tree/v3.4.0).

## Prerequisites
Download the required input files from the official
[ExaWind Benchmarks repository](https://github.com/Exawind/exawind-benchmarks/tree/main/amr-wind/atmospheric_boundary_layer/neutral/input_files) and place them in a folder named `SimulationFiles`.

Make sure to update the simulation configuration file by
changing the `time.stop_time` parameter to:

```
time.stop_time = 120.0
```

You're now ready to launch your simulation!

## Running the Simulation
Here's how to launch your AMR-Wind simulation on a 2-machine MPI cluster
using the Inductiva Python API:

```python
"""AMR-Wind example."""
import inductiva

# Allocate a multi-machine MPI cluster on Google Cloud Platform
cloud_machine = inductiva.resources.MPICluster(
    machine_type="c2d-highcpu-112",
    num_machines=2,
    data_disk_gb=200,
    spot=True
)

# Initialize the Simulator
amr_wind = inductiva.simulators.AmrWind(
    version="3.4.1"
)

# Run the simulation
task = amr_wind.run(
    input_dir="/Path/to/SimulationFiles",
    sim_config_filename="abl_neutral.inp",
    on=cloud_machine
)

# Wait for the simulation to finish and download results
task.wait()
cloud_machine.terminate()
task.download_outputs()
```

This setup lets you seamlessly scale your simulations across multiple nodes, delivering faster results and 
enabling you to tackle larger, more detailed CFD problems than ever before.

> For an example of running your simulation on a single machine, check out this [tutorial](quick-start).

## Results
Here are the results of running this simulation on a multi-node MPI cluster compared to the 
baseline single-node configuration.

| Machine Type    | Nº Machines | vCPUs | Duration (s) | Speedup |
| --------------- | ------------ | ----- | ------------ | ------- |
| c2d-highcpu-112 | 1            | 112   | 1472         | 1.00x   |
| c2d-highcpu-112 | 2            | 224   | 1404         | 1.05x   |
| c2d-highcpu-112 | 4            | 448   | 1461         | 1.01x   |
| c2d-highcpu-112 | 8            | 896   | 977          | 1.51x   |
| c2d-highcpu-112 | 12           | 1344  | 1408         | 1.05x   |
| c2d-highcpu-112 | 16           | 1792  | 1001         | 1.47x   |

For this simulation setup, the **fastest configuration** was achieved using a
cluster of **8 machines (896 vCPUs)**, completing the run in **977 seconds**.

Interestingly, increasing the number of machines beyond this point did
not improve performance. This is likely due to the diminishing
returns of parallelization for smaller or moderately sized problems, where the 
communication overhead between nodes begins to offset the gains from additional
compute resources.

To explore these results in more detail, check out our [benchmarks page](https://inductiva.ai/guides/amr-wind/mpi-cluster-benchmarks).


