# Run Your First Simulation
This tutorial will show you how to run SNL-SWAN simulations using the Inductiva API. 

We will cover the tutorial example from the test files folder of the [official SNL-SWAN documentation](https://sandialabs.github.io/SNL-SWAN/tutorial.html) to help you get started with simulations.

## Prerequisites
1. Download the required files [here](https://sandialabs.github.io/SNL-SWAN/_downloads/ExampleFiles.zip) and rename the `INPUT` file to `INPUT.swn`. 
2. To resolve [this issue](https://github.com/sandialabs/SNL-SWAN/issues/8), remove the following lines of code from the `INPUT.swn` file:
`OBSTACLE TRANS 0.3 REFL 0.00 LINE 400 400 400 450`
`OBSTACLE TRANS 0.3 REFL 0.00 LINE  450 500 450 550`

Once these steps are complete, you'll be ready to send your simulation to the Cloud.

## Running an SNL-SWAN Simulation
Here is the code required to run an SNL-SWAN simulation using the Inductiva API:

```python
"""SNL-SWAN example."""
import inductiva

# Allocate cloud machine on Google Cloud Platform
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c2d-highcpu-4",
	spot=True)

# Initialize the Simulator
snl_swan = inductiva.simulators.SNLSWAN( \
    version="2.2")

# Run simulation
task = snl_swan.run(
    input_dir="/Path/to/ExampleFiles",
    sim_config_filename="INPUT.swn",
    on=cloud_machine)

# Wait for the simulation to finish and download the results
task.wait()
cloud_machine.terminate()

task.download_outputs()

task.print_summary()
```

In this basic example, we're using a cloud machine (`c2d-highcpu-4`) equipped with 4 virtual CPUs. 
For larger or more compute-intensive simulations, consider adjusting the `machine_type` parameter to select 
a machine with more virtual CPUs and increased memory capacity. You can explore the full range of available machines [here](https://console.inductiva.ai/machine-groups/instance-types).

> **Note**: Setting `spot=True` enables the use of spot machines, which are available at substantial discounts. 
> However, your simulation may be interrupted if the cloud provider reclaims the machine.

To adapt this script for other SNL-SWAN simulations, replace `input_dir` with the
path to your SNL-SWAN input files and set the `sim_config_filename` accordingly.

When the simulation is complete, we terminate the machine, download the results and print a summary of the simulation as shown below.

```
Task status: Success

Timeline:
	Waiting for Input         at 22/04, 15:00:31      0.807 s
	In Queue                  at 22/04, 15:00:32      43.516 s
	Preparing to Compute      at 22/04, 15:01:15      2.381 s
	In Progress               at 22/04, 15:01:18      36.52 s
		├> 1.149 s         dd if=/dev/stdin of=machinefile
		└> 35.215 s        swanrun -input INPUT.swn -mpi 4
	Finalizing                at 22/04, 15:01:54      0.44 s
	Success                   at 22/04, 15:01:55      

Data:
	Size of zipped output:    71.71 KB
	Size of unzipped output:  892.21 KB
	Number of output files:   13

Estimated computation cost (US$): 0.00036 US$
```

As you can see in the "In Progress" line, the part of the timeline that represents the actual execution of the simulation, 
the core computation time of this simulation was approximately 36.5 seconds.

It's that simple!
