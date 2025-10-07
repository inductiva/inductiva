# Run Your First Simulation
This tutorial will show you how to run CalculiX simulations using the Inductiva API. 

We will cover the `large stuctural test example` from the [official CalculiX website](https://www.dhondt.de/) to help you get started with simulations.

## Prerequisites
Download the required files [here](https://www.dhondt.de/ccx_2.22.structest.tar.bz2) and the simulation files will be placed inside the `CalculiX/ccx_2.22/structest` folder. Then, you’ll be ready to send your simulation to the Cloud.

## Running a CalculiX Simulation
Here is the code required to run a CalculiX simulation using the Inductiva API:

```python
"""CalculiX example"""
import inductiva

# Instantiate machine group
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c2d-highcpu-4",
	spot=True)

# Initialize the Simulator
calculix = inductiva.simulators.CalculiX( \
    version="2.22")

# Run simulation with config files in the input directory
task = calculix.run( \
    input_dir="/Path/to/structest",
    sim_config_filename="contact2e.inp",
    on=cloud_machine)

task.wait()
cloud_machine.terminate()

task.download_outputs()

task.print_summary()

```

In this basic example, we're using a cloud machine (`c2d-highcpu-4`) equipped with 4 virtual CPUs. 
For larger or more compute-intensive simulations, consider adjusting the `machine_type` parameter to select 
a machine with more virtual CPUs and increased memory capacity. You can explore the full range of available machines [here](https://console.inductiva.ai/machine-groups/instance-types).

> **Note**: Setting `spot=True` enables the use of [spot machines](../how-it-works/machines/spot-machines.md), which are available at substantial discounts. 
> However, your simulation may be interrupted if the cloud provider reclaims the machine.

To adapt this script for other CalculiX simulations, replace `input_dir` with the
path to your CalculiX input files and set the `sim_config_filename` accordingly.

When the simulation is complete, we terminate the machine, download the results and print a summary of the simulation as shown below.

```
Task status: Success

Timeline:
	Waiting for Input         at 12/08, 11:29:31      0.86 s
	In Queue                  at 12/08, 11:29:32      42.632 s
	Preparing to Compute      at 12/08, 11:30:15      1.404 s
	In Progress               at 12/08, 11:30:16      353.357 s
		└> 353.201 s       ccx -i contact2e
	Finalizing                at 12/08, 11:36:09      0.514 s
	Success                   at 12/08, 11:36:10      

Data:
	Size of zipped output:    1.41 MB
	Size of unzipped output:  4.77 MB
	Number of output files:   9

Estimated Task Compute Cost = 0.0023 US$
Task Orchestration Fee = 0.01 US$
Total Estimated Cost = 0.0123 US$
Learn more about costs at: https://inductiva.ai/guides/how-it-works/basics/how-much-does-it-cost
```

As you can see in the "In Progress" line, the part of the timeline that represents the actual execution of the simulation, 
the core computation time of this simulation was 353 seconds (approximately 5 minutes and 53 seconds).

```{banner_small}
:origin: calculix_quick_start
```