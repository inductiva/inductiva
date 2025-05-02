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
    machine_type="c2d-highcpu-4",
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
	Waiting for Input         at 21/04, 16:02:53      1.589 s
	In Queue                  at 21/04, 16:02:55      34.667 s
	Preparing to Compute      at 21/04, 16:03:29      4.888 s
	In Progress               at 21/04, 16:03:34      27.232 s
		└> 27.087 s        /opt/openmpi/4.1.6/bin/mpirun --use-hwthread-cpus OpenSeesMP Example.tcl
	Finalizing                at 21/04, 16:04:01      1.139 s
	Success                   at 21/04, 16:04:03      

Data:
	Size of zipped output:    19.60 MB
	Size of unzipped output:  49.79 MB
	Number of output files:   120

Estimated computation cost (US$): 0.00031 US$
```

As you can see in the "In Progress" line, the part of the timeline that represents the actual execution of the simulation, 
the core computation time of this simulation was approximately 27.2 seconds.

Although it's short, there's still room for improvement to reduce the processing
time.

### Scaling Up Your Simulation  
Scaling up your simulation is as simple as updating the `machine_type` parameter to a 16 vCPU machine (`c2d-highcpu-16`).

By increasing the number of vCPUs, we've reduced the processing time from 27.2 
to 9.4 seconds.

Here are the results of running the same simulation on a few machines:

|  Machine Type  | Virtual CPUs |     Time     | Estimated Cost |
|:--------------:|:------------:|:------------:|:--------------:|
|  c2d-highcpu-4 |       4      | 27.2 seconds | 0.00031 US$    |
|  c2d-highcpu-8 |       8      | 14.2 seconds | 0.00034 US$    |
| c2d-highcpu-16 |      16      | 9.4 seconds  | 0.00046 US$    |

Still in the testing phase? No problem! Just skip this step for now and start
with a machine with fewer vCPUs. Once you're satisfied with your results, you
can seamlessly scale your OpenSees simulation.

## Run OpeenSeesPy
To run OpenSees scripts written in Python, all you need to do is change the `interface` parameter to `python` to match the
file type of your OpenSeesPy use case. 

It's that simple!
