# Run the Case on Inductiva

## Running the Simulation
The following Python script runs the simulation for the **2021 Spring deployment**. Save the following code as a `.py` file before running it.

```python
"""Delft3D Simulation."""
import inductiva

# Allocate a machine on Google Cloud Platform
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c2d-highcpu-32",
    spot=True)

# Initialize the Simulator
delft3d = inductiva.simulators.Delft3D(\
    version="6.04.00")

# Run simulation
task = delft3d.run( \
    input_dir="/Path/to/Spring2021_Delft3D_setup",
    commands = ["mpirun -np 32 d_hydro.exe config_d_hydro.xml"],
    on=cloud_machine)

# Wait for the simulation to finish and download the results
task.wait()
cloud_machine.terminate()

task.download_outputs()

task.print_summary()
```

In this example, we use a `c2d-highcpu-32` machine with 32 virtual CPUs. Its performance is slightly better than a high-end desktop. You can explore the full range of available machines [here](https://console.inductiva.ai/machine-groups/instance-types).

> **Note**: Setting `spot=True` enables the use of [spot machines](../how-it-works/machines/spot-machines.md), which are available at substantial discounts.
> However, your simulation may be interrupted if the cloud provider reclaims the machine.

When the simulation is complete, we terminate the machine, download the results and print a summary of the simulation as shown below.

```
Task status: Success

Timeline:
	Waiting for Input         at 22/10, 09:15:47      1.742 s
	In Queue                  at 22/10, 09:15:49      41.202 s
	Preparing to Compute      at 22/10, 09:16:30      3.284 s
	In Progress               at 22/10, 09:16:33      64214.631 s
		└> 64214.365 s     mpirun -np 32 d_hydro.exe config_d_hydro.xml
	Finalizing                at 23/10, 03:06:48      214.477 s
	Success                   at 23/10, 03:10:22

Data:
	Size of zipped output:    5.68 GB
	Size of unzipped output:  18.34 GB
	Number of output files:   64

Total estimated cost (US$): 2.79 US$
	Estimated computation cost (US$): 2.78 US$
	Task orchestration fee (US$): 0.010 US$

Note: A per-run orchestration fee (0.010 US$) applies to tasks run from 01 Dec 2025, in addition to the computation costs.
Learn more about costs at: https://inductiva.ai/guides/how-it-works/basics/how-much-does-it-cost
```

As you can see in the "In Progress" line, the part of the timeline that represents the actual execution of the simulation, the core computation time of this simulation was approximately **hours**.

To run the 2021 Fall deployment, simply update the `input_dir` parameter accordingly.

## Scaling Up Your Simulation
Scaling up is simple: either increase the number of vCPUs by choosing a larger `machine_type` (e.g., from `c2d-highcpu-32` to `c2d-highcpu-56`) or switch to a machine from the latest-generation **c4d series**. Both approaches can significantly reduce runtime.

⚠️ **Important**: If you change the number of vCPUs, be sure to update the `mpirun` command accordingly (e.g., `-np 56`) to match your machine configuration.

Below are the results of running the **2021 Spring deployment** across different computational setups:

| Machine Type      | vCPUs | Execution Time | Estimated Cost (USD)|
|-------------------|-------|----------------|---------------------|
| c2d-highcpu-32    | 32    | 17h, 54 min    | 2.78                |
| c2d-highcpu-56    | 56    | 13h, 14 min    | 3.52                |
| c2d-highcpu-112   | 112   | 10h, 33 min    | 5.52                |
| c4d-highcpu-48    | 48    | 9h, 30 min     | 7.93                |
| c4d-highcpu-96    | 96    | 6h, 49 min     | 11.32               |

As shown above, increasing the number of vCPUs or switching to newer-generation machines leads to a clear reduction in runtime, though at a higher computational cost. The optimal setup depends on your priorities, whether minimizing time-to-results or optimizing cost-efficiency. The **Inductiva** platform makes it easy to explore this trade-off, enabling you to scale Delft3D simulations seamlessly across different compute configurations.

```{banner_small}
:origin: swan
```
