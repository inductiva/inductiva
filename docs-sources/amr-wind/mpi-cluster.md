# Running AMR-Wind Simulations Across Multiple Machines Using MPI

As computational fluid dynamics (CFD) problems grow in complexity and scale, so
does the demand for computing power. **AMR-Wind** was built to scale efficiently
across thousands of cores. However, even the most powerful single machine can
eventually become a bottleneck, particularly when dealing with large domains or
fine mesh resolutions.

To overcome these limitations, you can run AMR-Wind simulations across
**multiple machines** using an **MPI (Message Passing Interface) cluster**. This
enables you to combine the CPU resources of several machines, significantly
increasing the total number of virtual CPUs (vCPUs) available to your simulation
and accelerating runtime.

With **Inductiva**, running your simulations on a multi-node MPI cluster becomes
as easy as using any other single node option.

**Changing MachineGroup to MPICluster**
```diff
# Allocate a multi-machine MPI cluster on Google Cloud Platform
- cloud_machine = inductiva.resources.MachineGroup(
+ cloud_machine = inductiva.resources.MPICluster(
    machine_type="c2d-highcpu-112",
+   num_machines=2,
    data_disk_gb=100,
    spot=True)
```

## Prerequisites

To begin, download the necessary input files from the official
[Exawind Benchmarks repository](https://github.com/Exawind/exawind-benchmarks/tree/main/amr-wind/atmospheric_boundary_layer/neutral/input_files) and place them in a folder named `SimulationFiles`.

Before running, make sure to update the simulation configuration file by
changing the `time.stop_time` parameter to:

```
time.stop_time = 120.0
```

You're now ready to launch your simulation.

---

## Running the Simulation

Here's how to launch your AMR-Wind simulation on a 2 machines MPI cluster
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

# Initialize the simulator
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

This setup allows you to seamlessly scale your simulations across multiple
nodes, getting results faster, and tackling larger and more detailed CFD problems
than ever before.

You can find an example of how to run your simulation on a single machine [here](quick-start).

## Results

Here are the results of running this simulation on multiple MPI cluster compared
with the baseline of running in a single node configuration.

| Machine Type    | Num Machines | vCPUs | Duration (s) | Speedup |
| --------------- | ------------ | ----- | ------------ | ------- |
| c2d-highcpu-112 | 1            | 112   | 1472         | 1.00x   |
| c2d-highcpu-112 | 2            | 224   | 1404         | 1.05x   |
| c2d-highcpu-112 | 4            | 448   | 1461         | 1.01x   |
| c2d-highcpu-112 | 8            | 896   | 977          | 1.51x   |
| c2d-highcpu-112 | 12           | 1344  | 1408         | 1.05x   |
| c2d-highcpu-112 | 16           | 1792  | 1001         | 1.47x   |


For this simulation setup, the **fastest configuration** was achieved using a
cluster with **8 machines (896 vCPUs)**, completing in **977 seconds**.
Interestingly, increasing the number of machines beyond this point did
**not lead to better performance**. This can be attributed to the diminishing
returns of parallelization for smaller or moderately sized problems, where
communication overhead between nodes begins to offset the gains from additional
compute resources.

## What Happens When We Increase the Simulation Complexity?

To investigate how performance scales with problem size, we increased the
simulation complexity by refining the grid resolution. Specifically, we changed:

```
amr.n_cell              = 512 512 184 # Grid cells at coarsest AMRlevel
```

to

```
amr.n_cell              = 1024 1024 368 # Grid cells at coarsest AMRlevel
```

This effectively **increased the total number of grid cells by a factor of 8**,
leading to a significantly more computationally demanding simulation. In these
cases, **larger clusters with more vCPUs are expected to perform better**, as
the computational workload can be more effectively distributed across machines.
The communication overhead becomes relatively smaller compared to the total
compute time, which improves scalability.

We'll now compare the performance across the same cluster configurations under
this higher-resolution setup.

| Machine Type    | Num Machines | Cores | Duration (s) | Speedup |
| --------------- | ------------ | ----- | ------------ | ------- |
| c2d-highmem-112 | 1            | 112   | 8640         | 1.00x   |
| c2d-highmem-112 | 2            | 224   | 5400         | 1.60x   |
| c2d-highmem-112 | 4            | 448   | 3550         | 2.43x   |
| c2d-highmem-112 | 8            | 896   | 2532         | 3.41x   |
| c2d-highmem-112 | 12           | 1344  | 2055         | 4.20x   |
| c2d-highmem-112 | 16           | 1792  | 1848         | 4.68x   |


The results show that, with the higher-resolution grid, simulation runtimes
decreased consistently as more machines were added, unlike in the lower-resolution
case, where gains plateaued beyond 8 machines. The best performance was observed
using 16 machines (1792 vCPUs), bringing the runtime down to just 1848 seconds,
a significant improvement compared to the single-machine setup at 8640 seconds.

However, the speedup does not scale linearly with the number of vCPUs. While the
4 and 8 machine setups provided solid performance gains, adding more machines
beyond that point yielded diminishing returns.
