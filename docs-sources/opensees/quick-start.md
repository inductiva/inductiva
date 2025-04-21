# Run Your First Simulation
This tutorial will show you how to run OpenSees simulations using the Inductiva API. 

This tutorial will cover the `SmallMP` use case of the OpenSees Examples, available in the [OpenSees GitHub repository](https://github.com/OpenSees/OpenSees), to help you get started with simulations.

We will also demonstrate Inductiva’s ability to efficiently scale this use case on a more powerful machine.

## Prerequisites
Download the required files [here](https://github.com/OpenSees/OpenSees/tree/master/EXAMPLES/SmallMP) and place them in a folder called `SmallMP`. Then, you’ll be ready to send your simulation to the Cloud.

## Running an OpenSees Simulation
Here is the code required to run OpenSees simulation using the Inductiva API:

```python
"""OpenSees Simulation."""
import inductiva

# Allocate a machine on Google Cloud Platform
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c2-standard-4",
	spot=True)

# Initialize the Simulator
opensees = inductiva.simulators.OpenSees( \
    interface="tcl",
    version="3.7.1")

# Run simulation
task = opensees.run( \
    input_dir="/Path/to/SmallMP",
    sim_config_filename="Example.tcl",
    on=cloud_machine)

# Wait for the simulation to finish and download the results
task.wait()
cloud_machine.terminate()

task.download_outputs()

task.print_summary()
```

To adapt this script for other OpenSees simulations, replace `input_dir` with the
path to your OpenSees input files and set the `sim_config_filename` accordingly.

When the simulation is complete, we terminate the machine, download the results and print a summary of the simulation as shown below.

```
Task status: Success

Timeline:
	Waiting for Input         at 25/02, 16:35:31      1.995 s
	In Queue                  at 25/02, 16:35:33      20.562 s
	Preparing to Compute      at 25/02, 16:35:54      3.996 s
	In Progress               at 25/02, 16:35:58      28.221 s
		└> 28.09 s         /opt/openmpi/4.1.6/bin/mpirun --use-hwthread-cpus --np 4 OpenSeesMP Example.tcl
	Finalizing                at 25/02, 16:36:26      1.339 s
	Success                   at 25/02, 16:36:27      

Data:
	Size of zipped output:    23.23 MB
	Size of unzipped output:  59.07 MB
	Number of output files:   186

Estimated computation cost (US$): 0.00066 US$
```

As you can see in the "In Progress" line, the part of the timeline that represents the actual execution of the simulation, 
the core computation time of this simulation was 28.2 seconds.

Although it's short, there's still room for improvement to reduce the processing
time.

### Scaling Up Your Simulation  

Scaling up your simulation is as simple as updating the `machine_type` parameter to a 16 vCPU machine (`c2-standard-16`).

By increasing the number of vCPUs, we've reduced the processing time from 28.2 
to 10.2 seconds.

Here are the results of running the same simulation on a few machines:

|  Machine Type  | Virtual CPUs |     Time     | Estimated Cost |
|:--------------:|:------------:|:------------:|:--------------:|
|  c2-standard-4 |       4      | 28.2 seconds | 0.00066 US$    |
|  c2-standard-8 |       8      | 15.2 seconds | 0.00094 US$    |
| c2-standard-16 |      16      | 10.2 seconds | 0.0011 US$     |

Still in the testing phase? No problem! Just skip this step for now and start
with a machine with fewer vCPUs. Once you're satisfied with your results, you
can seamlessly scale your OpenSees simulation.

## Run OpeenSeesPy
To run OpenSees scripts written in Python, all you need to do is change the `interface` parameter to `python` to match the
file type of your OpenSeesPy use case. 

It's that simple!
