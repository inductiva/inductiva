# Run Your First Simulation
This tutorial will show you how to run SNL-SWAN simulations using the Inductiva API. 

We will cover the tutorial example from the test files folder of the [official SNL-SWAN documentation](https://sandialabs.github.io/SNL-SWAN/tutorial.html) to help you get started with simulations.

## Prerequisites
Download the required files [here](https://sandialabs.github.io/SNL-SWAN/_downloads/ExampleFiles.zip). Then, you’ll be ready to send your simulation to the Cloud.


## Running an SNL-SWAN Simulation
Here is the code required to run an SNL-SWAN simulation using the Inductiva API:

```python
"""SNL-SWAN example."""
import inductiva

# Allocate cloud machine on Google Cloud Platform
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c2d-standard-4",
	spot=True)

# Initialize the Simulator
snl_swan = inductiva.simulators.SNLSWAN( \
    version="2.2")

# Run simulation
task = snl_swan.run(
    input_dir="/Path/to/ExampleFiles",
    sim_config_filename="INPUT",
    on=cloud_machine)

# Wait for the simulation to finish and download the results
task.wait()
cloud_machine.terminate()

task.download_outputs()

task.print_summary()
```

> **Note**: `spot` machines are a lot cheaper but may be terminated by the provider if necessary.

To adapt this script for other SNL-SWAN simulations, replace `input_dir` with the
path to your SNL-SWAN input files and set the `sim_config_filename` accordingly.

When the simulation is complete, we terminate the machine, download the results and print a summary of the simulation as shown below.

```
Task status: Success

Timeline:
	Waiting for Input         at 07/04, 14:31:10      0.824 s
	In Queue                  at 07/04, 14:31:11      17.849 s
	Preparing to Compute      at 07/04, 14:31:28      2.056 s
	In Progress               at 07/04, 14:31:30      4.337 s
		├> 1.108 s         dd if=/dev/stdin of=machinefile
		└> 3.059 s         swanrun -input INPUT
	Finalizing                at 07/04, 14:31:35      0.921 s
	Success                   at 07/04, 14:31:36      

Data:
	Size of zipped output:    3.04 KB
	Size of unzipped output:  11.84 KB
	Number of output files:   6

Estimated computation cost (US$): 0.000085 US$
```

As you can see in the "In Progress" line, the part of the timeline that represents the actual execution of the simulation, 
the core computation time of this simulation was 4.3 seconds.

It's that simple!
