# Run Your First Simulation
This tutorial will show you how to run SWAN simulations using the Inductiva API. 

We will cover the `ring` use case from the [official SWAN documentation](https://swanmodel.sourceforge.io/download/download.htm) to help you get started with simulations.

## Prerequisites
Download the required files [here](https://swanmodel.sourceforge.io/download/zip/ring.tar.gz) and place them in a folder called `Ring`. Then, you’ll be ready to send your simulation to the Cloud.

## Running an SWAN Simulation
Here is the code required to run a SWAN simulation using the Inductiva API:

```python
"""SWAN example"""
import inductiva

# Allocate cloud machine on Google Cloud Platform
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c2d-highcpu-4",
	spot=True)

# Initialize the Simulator
swan = inductiva.simulators.SWAN(\
    version="41.51")

# Run simulation
task = swan.run(input_dir="/Path/to/Ring",
    sim_config_filename="ring.swn",
    on=cloud_machine)

# Wait for the simulation to finish and download the results
task.wait()
cloud_machine.terminate()

task.download_outputs()

task.print_summary()
```

> **Note**: `spot` machines are a lot cheaper but may be terminated by the provider if necessary.

To adapt this script for other SWAN simulations, replace `input_dir` with the
path to your SWAN input files and set the `sim_config_filename` accordingly.

When the simulation is complete, we terminate the machine, download the results and print a summary of the simulation as shown below.

```
Task status: Success

Timeline:
	Waiting for Input         at 21/04, 18:56:15      1.571 s
	In Queue                  at 21/04, 18:56:17      38.79 s
	Preparing to Compute      at 21/04, 18:56:56      2.99 s
	In Progress               at 21/04, 18:56:59      55.349 s
		├> 1.095 s         dd if=/dev/stdin of=machinefile
		└> 54.101 s        swanrun -input ring.swn -mpi 4
	Finalizing                at 21/04, 18:57:54      0.416 s
	Success                   at 21/04, 18:57:54      

Data:
	Size of zipped output:    118.60 KB
	Size of unzipped output:  488.62 KB
	Number of output files:   19

Estimated computation cost (US$): 0.00053 US$
```

As you can see in the "In Progress" line, the part of the timeline that represents the actual execution of the simulation, 
the core computation time of this simulation was 55.3 seconds.

It's that simple!
