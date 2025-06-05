# Run Your First Simulation
This tutorial will show you how to run Deft3D simulations using the Inductiva API. 

We will cover the `01_standard` use case from the examples available in the [Delft3D Subversion repository](https://svn.oss.deltares.nl/repos/delft3d/branches/releases/7545/),

## Prerequisites
Download the required files [here](https://svn.oss.deltares.nl/repos/delft3d/branches/releases/7545/examples/01_standard/) and save them to a folder named `SimulationFiles`. Then, you’ll be ready to
send your simulation to the Cloud.

## Running a Delft3D Simulation
Here is the code required to run the Delft3D simulation using the Inductiva API:

```python
"""Delft3D Simulation."""
import inductiva

# Allocate a machine on Google Cloud Platform
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c2-standard-4",
	spot=True)

# Initialize the Simulator
delft3d = inductiva.simulators.Delft3D(\
    version="6.04.00")

# Run simulation
task = delft3d.run( \
    input_dir="/Path/to/SimulationFiles",
    commands = ["mpirun -np 4 d_hydro.exe config_d_hydro.xml"],
    on=cloud_machine)

# Wait for the simulation to finish and download the results
task.wait()
cloud_machine.terminate()

task.download_outputs()

task.print_summary()
```

> **Note**: `spot` machines are available at substantial discounts, but your simulation job may be preempted if
> the Cloud provider reclaims the spot machine.

To adapt this script for other Delft3D simulations, replace `input_dir` with the
path to your Delft3D input files and set the the `commands` accordingly.

When the simulation is complete, we terminate the machine, download the results and print a summary of the simulation as shown below.

```
Task status: Success

Timeline:
	Waiting for Input         at 10/04, 14:17:59      0.695 s
	In Queue                  at 10/04, 14:18:00      36.277 s
	Preparing to Compute      at 10/04, 14:18:36      2.832 s
	In Progress               at 10/04, 14:18:39      5.191 s
		└> 5.073 s         mpirun -np 4 d_hydro.exe config_d_hydro.xml
	Finalizing                at 10/04, 14:18:44      0.449 s
	Success                   at 10/04, 14:18:45      

Data:
	Size of zipped output:    433.45 KB
	Size of unzipped output:  1.04 MB
	Number of output files:   15

Estimated computation cost (US$): 0.00017 US$

Go to https://console.inductiva.ai/tasks/iidkjkpk77yb79cq4w5qpgdq8 for more details.
```

As you can see in the "In Progress" line, the part of the timeline that represents the actual execution of the simulation, 
the core computation time of this simulation was approximately 5.2 seconds.

It's that simple!
