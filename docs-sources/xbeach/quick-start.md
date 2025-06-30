# Run Your First Simulation
This tutorial will show you how to run XBeach simulations using the Inductiva API. 

We will cover the `Bijleveld_surfbeat_s200` test case from the [official XBeach Subversion repository](https://svn.oss.deltares.nl/repos/xbeach/) to help you get started with simulations.

## Prerequisites
Download all of the required files [here](https://svn.oss.deltares.nl/repos/xbeach/testcases/Wong2016/Bijleveld_surfbeat_s200/) and place them into a folder named `SimulationFiles`. Then, you’ll be ready to send your simulation to the Cloud.

## Adjust Simulation Parameters (Optional)
For the purposes of this tutorial, we'll shorten the simulation time by a factor of 10. To do this, modify the `tstop` parameter in the `params.txt` file to 1800.

## Running an XBeach Simulation
Here is the code required to run an XBeach simulation using the Inductiva API:

```python
"""XBeach example"""
import inductiva

# Allocate cloud machine on Google Cloud Platform
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c2d-highcpu-16",
	spot=True)

# Initialize the Simulator
xbeach = inductiva.simulators.XBeach( \
    version="1.24")

# Run simulation
task = xbeach.run(input_dir="/Path/to/SimulationFiles",
    sim_config_filename="params.txt",
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

> **Note**: Setting `spot=True` enables the use of spot machines, which are available at substantial discounts. 
> However, your simulation may be interrupted if the cloud provider reclaims the machine.

To adapt the code for this or any other use case, simply replace `input_dir` with the path to your XBeach files before executing it in a Python script.

When the simulation is complete, we terminate the machine, download the results and print a summary of the simulation as shown below.

```
Task status: Success

Timeline:
	Waiting for Input         at 21/04, 19:24:34      1.145 s
	In Queue                  at 21/04, 19:24:35      38.265 s
	Preparing to Compute      at 21/04, 19:25:13      1.226 s
	In Progress               at 21/04, 19:25:14      368.298 s
		└> 368.152 s       /opt/openmpi/4.1.6/bin/mpirun --use-hwthread-cpus xbeach params.txt
	Finalizing                at 21/04, 19:31:22      5.112 s
	Success                   at 21/04, 19:31:28      

Data:
	Size of zipped output:    132.13 MB
	Size of unzipped output:  175.42 MB
	Number of output files:   13

Estimated computation cost (US$): 0.012 US$
```

As you can see in the "In Progress" line, the part of the timeline that represents the actual execution of the simulation, 
the core computation time of this simulation was 368.3 seconds (approximately 6 minutes and 8 seconds).

It's that simple!