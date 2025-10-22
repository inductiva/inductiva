# Run the Case on Inductiva

## Running the Baseline S0 Scenario
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
<Add>

Estimated Task Compute Cost =  US$
Task Orchestration Fee = 0.01 US$
Total Estimated Cost = US$
Learn more about costs at: https://inductiva.ai/guides/how-it-works/basics/how-much-does-it-cost
```

As you can see in the "In Progress" line, the part of the timeline that represents the actual execution of the simulation, the core computation time of this simulation was approximately **hours**.

To run the 2021 Fall deployment, simply update the `input_dir` parameter accordingly.

## Scaling Up Your Simulation  
Scaling up is simple — either increase the number of vCPUs by choosing a larger `machine_type` (e.g., from `c2d-highcpu-32` to `c2d-highcpu-56`) or switch to a machine from the latest-generation **c4d series**. Both approaches can significantly reduce runtime.

⚠️ **Important**: If you change the number of vCPUs, be sure to update the `mpirun` command accordingly (e.g., `-np 56`) to match your machine configuration.

| Scenario | Machine Type      | vCPUs | Execution Time | Estimated Cost (US$) |
|----------|-------------------|-------|----------------|---------------------|

<Conclusion>

```{banner_small}
:origin: swan
```