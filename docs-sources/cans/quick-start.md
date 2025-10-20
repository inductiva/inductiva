# Run Your First Simulation
This tutorial will show you how to run CaNS simulations using the Inductiva API. 

We will cover the `closed_box` use case from the examples available in the [official CaNS documentation](https://github.com/CaNS-World/CaNS), to help you get started with simulations.

## Prerequisites
1. Download the required files [here](https://github.com/CaNS-World/CaNS/tree/v2.4.0/examples/closed_box) and save them to a folder named `SimulationFiles`.
2. In the `SimulationFiles` folder, create a new folder named `data`. The simulator writes files to this folder and will encounter an error if it is not present.

Then, you’ll be ready to send your simulation to the Cloud.

## Running a CaNS Simulation
Here is the code required to run a CaNS simulation using the Inductiva API:

```python
"""CaNS example."""
import inductiva

# Allocate cloud machine on Google Cloud Platform
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c2d-highcpu-16",
	spot=True)

# Initialize the Simulator
cans = inductiva.simulators.CaNS(\
    version="2.4.0")

# Run simulation
task = cans.run(input_dir="/Path/to/SimulationFiles",
    sim_config_filename="input.nml",
    on=cloud_machine)

# Wait for the simulation to finish and download the results
task.wait()
cloud_machine.terminate()

task.download_outputs()

task.print_summary()
```

In this basic example, we're using a cloud machine (`c2d-highcpu-16`) equipped with 16 virtual CPUs. 
For larger or more compute-intensive simulations, consider adjusting the `machine_type` parameter to select 
a machine with more virtual CPUs and increased memory capacity. You can explore the full range of available machines [here](https://console.inductiva.ai/machine-groups/instance-types).

> **Note**: Setting `spot=True` enables the use of [spot machines](../how-it-works/machines/spot-machines.md), which are available at substantial discounts. 
> However, your simulation may be interrupted if the cloud provider reclaims the machine.

To adapt the code for this or any other use case, simply replace `input_dir` with the path to your CaNS input files and 
set the `sim_config_filename` accordingly. Be sure to specify the CaNS version compatible with your input files.

When the simulation is complete, we terminate the machine, download the results and print a summary of the simulation as shown below.

```
Task status: Success

Timeline:
	Waiting for Input         at 17/04, 15:25:52      0.797 s
	In Queue                  at 17/04, 15:25:52      73.225 s
	Preparing to Compute      at 17/04, 15:27:06      1.518 s
	In Progress               at 17/04, 15:27:07      126.298 s
		└> 126.177 s       /opt/openmpi/4.1.6/bin/mpirun --use-hwthread-cpus cans input.nml
	Finalizing                at 17/04, 15:29:13      1.236 s
	Success                   at 17/04, 15:29:15      

Data:
	Size of zipped output:    3.50 MB
	Size of unzipped output:  203.68 MB
	Number of output files:   3260

Estimated Task Compute Cost = 0.0042 US$
Task Orchestration Fee = 0.01 US$
Total Estimated Cost = 0.0142 US$
Learn more about costs at: https://inductiva.ai/guides/how-it-works/basics/how-much-does-it-cost
```

As you can see in the "In Progress" line, the part of the timeline that represents the actual execution of the simulation, 
the core computation time of this simulation was approximately 126.3 seconds (around 2 minutes).

```{banner_small}
:origin: cans
```