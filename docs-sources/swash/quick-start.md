# Run Your First Simulation
This tutorial will show you how to run SWASH simulations using the Inductiva API. 

We will cover the `Berkhoff shoal` use case from the [official SWASH documentation](https://swash.sourceforge.io/), to help you get started with simulations.

## Prerequisites
Download the SWASH test case collection [here](https://swash.sourceforge.io/). You'll find the input file needed to run this use case in the `l41berkh` folder. Then, you'll be ready to send your simulation to the Cloud!

## Running a SWASH Simulation
Here is the code required to run a SWASH simulation using the Inductiva API:

```python
"""SWASH example."""
import inductiva

# Allocate cloud machine on Google Cloud Platform
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c3d-standard-16",
	spot=True)

# Initialize the Simulator
swash = inductiva.simulators.SWASH(\
    version="11.01")

# Run simulation
task = swash.run(input_dir="/Path/to/l41berkh",
    sim_config_filename="l41ber02.sws",
    on=cloud_machine)

# Wait for the simulation to finish and download the results
task.wait()
cloud_machine.terminate()

task.download_outputs()

task.print_summary()
```

> **Note**: `spot` machines are a lot cheaper but may be terminated by the provider if necessary.

To adapt the code for this or any other use case, simply replace `input_dir` with the path to your SWASH input files and 
set the `sim_config_filename` accordingly.

When the simulation is complete, we terminate the machine, download the results and print a summary of the simulation as shown below.

```
Task status: Success

Timeline:
	Waiting for Input         at 03/04, 11:34:37      0.928 s
	In Queue                  at 03/04, 11:34:38      17.582 s
	Preparing to Compute      at 03/04, 11:34:56      2.01 s
	In Progress               at 03/04, 11:34:58      42.406 s
		├> 1.139 s         dd if=/dev/stdin of=machinefile
		└> 41.12 s         swashrun -input l41ber01.sws -mpi 16
	Finalizing                at 03/04, 11:35:40      1.581 s
	Success                   at 03/04, 11:35:42      

Data:
	Size of zipped output:    35.68 MB
	Size of unzipped output:  113.87 MB
	Number of output files:   31

Estimated computation cost (US$): 0.0024 US$
```

As you can see in the "In Progress" line, the part of the timeline that represents the actual execution of the simulation, 
the core computation time of this simulation was 42.4 seconds.

It's that simple!
