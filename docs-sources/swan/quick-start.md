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
    machine_type="c3d-standard-16",
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
	Waiting for Input         at 08/04, 17:29:30      1.013 s
	In Queue                  at 08/04, 17:29:31      44.577 s
	Preparing to Compute      at 08/04, 17:30:16      5.274 s
	In Progress               at 08/04, 17:30:21      20.31 s
		├> 1.059 s         dd if=/dev/stdin of=machinefile
		└> 19.08 s         swanrun -input ring.swn -mpi 16
	Finalizing                at 08/04, 17:30:42      0.477 s
	Success                   at 08/04, 17:30:42      

Data:
	Size of zipped output:    175.63 KB
	Size of unzipped output:  639.49 KB
	Number of output files:   43

Estimated computation cost (US$): 0.0014 US$
```

As you can see in the "In Progress" line, the part of the timeline that represents the actual execution of the simulation, 
the core computation time of this simulation was 20.3 seconds.

It's that simple!
