# test
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
    machine_type="c2d-highcpu-4",
	spot=True)

# Initialize the Simulator
swash = inductiva.simulators.SWASH(\
    version="11.01")

# Run simulation
task = swash.run(input_dir="/Path/to/l41berkh",
    sim_config_filename="l41ber01.sws",
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
	Waiting for Input         at 21/04, 19:43:50      0.941 s
	In Queue                  at 21/04, 19:43:51      35.635 s
	Preparing to Compute      at 21/04, 19:44:27      2.26 s
	In Progress               at 21/04, 19:44:29      140.425 s
		├> 1.063 s         dd if=/dev/stdin of=machinefile
		└> 139.205 s       swashrun -input l41ber01.sws -mpi 4
	Finalizing                at 21/04, 19:46:49      1.665 s
	Success                   at 21/04, 19:46:51      

Data:
	Size of zipped output:    35.41 MB
	Size of unzipped output:  110.91 MB
	Number of output files:   19

Estimated computation cost (US$): 0.0013 US$
```

As you can see in the "In Progress" line, the part of the timeline that represents the actual execution of the simulation, 
the core computation time of this simulation was 140.4 seconds (approximately 2 minutes and 20 seconds).

It's that simple!
