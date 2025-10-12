# Running the Scenarios on Inductiva

## Running the Baseline S0 Scenario
The following Python script runs the simulation for the configuration `S0/polynya2D_20211007`.swn. Save this code into a `.py` file and place it inside the `tutorial` directory before running it.

```python
"""SWAN example"""
import inductiva

# Allocate cloud machine on Google Cloud Platform
cloud_machine = inductiva.resources.MachineGroup( \
	provider="GCP",
	machine_type="c4d-highcpu-32",
	data_disk_gb=20,
	spot=True)

# Initialize the Simulator
swan = inductiva.simulators.SWAN()

# Run simulation
task = swan.run(input_dir="/Path/to/polynya",
	sim_config_filename="S0/polynya2D_20211007.swn",
	on=cloud_machine)

# Wait for the simulation to finish and download the results
task.wait()
cloud_machine.terminate()

task.print_summary()
```

In this example, we use a `c4d-highcpu-32` machine, based on 5th Gen AMD EPYC processors (2024), equipped with 32 virtual CPUs. You can explore the full range of available machines [here](https://console.inductiva.ai/machine-groups/instance-types).

> **Note**: Setting `spot=True` enables the use of [spot machines](../how-it-works/machines/spot-machines.md), which are available at substantial discounts. 
> However, your simulation may be interrupted if the cloud provider reclaims the machine.

When the simulation is complete, we terminate the machine, download the results and print a summary of the simulation as shown below.

```
Task status: Success

Timeline:
	Waiting for Input         at 05/10, 10:54:03      1.262 s
	In Queue                  at 05/10, 10:54:04      50.123 s
	Preparing to Compute      at 05/10, 10:54:54      6.886 s
	In Progress               at 05/10, 10:55:01      12956.885 s
		├> 1.068 s         dd if=/dev/stdin of=machinefile
		└> 12955.593 s     swanrun -input polynya2D_20211007.swn -mpi 32
	Finalizing                at 05/10, 14:30:58      6.567 s
	Success                   at 05/10, 14:31:04      

Data:
	Size of zipped output:    930.15 MB
	Size of unzipped output:  1.92 GB
	Number of output files:   180

Estimated Task Compute Cost = 2.01 US$
Task Orchestration Fee = 0.01 US$
Total Estimated Cost = 2.02 US$
Learn more about costs at: https://inductiva.ai/guides/how-it-works/basics/how-much-does-it-cost
```

As you can see in the "In Progress" line, the part of the timeline that represents the actual execution of the simulation, the core computation time of this simulation was approximately **3.6 hours**.

To run the `S2_f5` or other scenarios, simply update the `sim_config_filename` parameter accordingly.

## Scaling Up Your Simulation  
Scaling up is straightforward — just increase the `machine_type` to a machine with more virtual CPUs. For example, doubling the vCPUs from `c4d-highcpu-32` to `c4d-highcpu-64` can significantly reduce runtime.

| Scenario | Machine Type      | vCPUs | Execution Time | Estimated Cost (US$) |
|----------|-------------------|-------|----------------|---------------------|
| **S0**   | c4d-highcpu-32    | 32    | 3h 36m         | 2.00                |
|          | c4d-highcpu-64    | 64    | 2h 26m         | 2.70                |
| **S2_f5**| c4d-highcpu-32    | 32    | 2h 37m         | 1.45                |
|          | c4d-highcpu-48    | 48    | 1h 58m         | 1.64                |
|          | c4d-highcpu-64    | 64    | 1h 42m         | 1.89                |
|          | c4d-highcpu-96    | 96    | 1h 34m         | 2.60                |

We observe a significant, though not strictly linear, improvement in execution time as the number of virtual CPUs increases, accompanied by a corresponding rise in cost. For larger simulations or tight deadlines, opting for more powerful machines can yield substantial time savings.

```{banner_small}
:origin: swan
```