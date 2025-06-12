# Run AMR-Wind Simulations Across Multiple Machines Using MPI
As computational fluid dynamics (CFD) problems increase in complexity and scale, so
does the demand for computing power. **AMR-Wind** is designed to scale efficiently across thousands of cores, but even the most powerful single machine can become a bottleneck, particularly when handling large domains or very fine mesh resolutions.

To overcome these limitations, you can run AMR-Wind simulations across
**multiple machines** using an **MPI (Message Passing Interface) cluster**. This approach combines the 
CPU resources of several machines, significantly boosting the total number of virtual CPUs (vCPUs) 
available and accelerating your simulation’s runtime.

With **Inductiva**, running your simulations on a multi-node MPI cluster is as straightforward as 
using any single-node option.

All you need to do is update the resource allocation from `MachineGroup` to `MPICluster`, as follows:

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

Make sure to update the simulation configuration file (`abl_neutral.inp`) by
setting the `time.stop_time` parameter to:

```
time.stop_time = 120.0
```

To evaluate scalability under heavier workloads, increase the mesh refinement by modifying the `amr.n_cell` parameter in the same file. Change it from:

```
amr.n_cell              = 512 512 184 # Grid cells at coarsest AMRlevel
```

to:

```
amr.n_cell              = 1024 1024 368 # Grid cells at coarsest AMRlevel
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

## Results
Below are the results of running this simulation on a multi-node MPI cluster, compared to the single-node configuration (shown in the first row):

<table>
  <tr>
    <td>Machine Type</td>
    <td>Nº of Machines</td>
    <td>Cores</td>
    <td>Duration (min:s)</td>
    <td>Speedup</td>
  </tr>
  <tr>
    <td>c2d-highmem-112</td>
    <td>1</td>
    <td>112</td>
    <td>144:00</td>
    <td>Baseline</td>
  </tr>
  <tr>
    <td>c2d-highmem-112</td>
    <td>2</td>
    <td>224</td>
    <td>90:00</td>
    <td>1.60x</td>
  </tr>
  <tr>
    <td>c2d-highmem-112</td>
    <td>4</td>
    <td>448</td>
    <td>59:10</td>
    <td>2.43x</td>
  </tr>
  <tr>
    <td>c2d-highmem-112</td>
    <td>8</td>
    <td>896</td>
    <td>42:12</td>
    <td>3.41x</td>
  </tr>
  <tr>
    <td>c2d-highmem-112</td>
    <td>12</td>
    <td>1344</td>
    <td>34:15</td>
    <td>4.20x</td>
  </tr>
  <tr>
    <td>c2d-highmem-112</td>
    <td>16</td>
    <td>1792</td>
    <td>30:48</td>
    <td>4.68x</td>
  </tr>
</table>

Runtime consistently decreased as the cluster size grew. The best performance was recorded with **16 machines** (1792 vCPUs), completing the simulation in **30 minutes and 48 seconds** — a major improvement over the single-machine runtime of **144 minutes**.

Nonetheless, speedup does not scale linearly with the number of vCPUs. While the 4- and 8- machine setups showed solid performance gains, adding more machines beyond that point yielded diminishing returns.