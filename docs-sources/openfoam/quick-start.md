# Run Your First Simulation
This tutorial will show you how to run OpenFOAM simulations using the Inductiva API. 

We will cover the `motorbike` use case from the [OpenFOAM Foundation GitHub repository](https://github.com/OpenFOAM/OpenFOAM-8/tree/version-8/tutorials), to help you get started with simulations.

## Prerequisites
Download the required files [here](https://github.com/OpenFOAM/OpenFOAM-8/tree/version-8/tutorials/incompressible/simpleFoam/motorBike) and place them in a folder called `SimulationFiles`. Then, you’ll be ready to send your simulation to the Cloud.

## Running an OpenFOAM Simulation
Here is the code required to run an OpenFOAM simulation using the Inductiva API:

```python
"""OpenFOAM example"""
import inductiva

# Allocate cloud machine on Google Cloud Platform
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c2d-highcpu-16",
	spot=True)

# Initialize the Simulator
OpenFOAM = inductiva.simulators.OpenFOAM( \
    version="8",
	distribution="foundation")

# Run simulation
task = OpenFOAM.run(input_dir="/Path/to/SimulationFiles",
    shell_script="./Allrun",
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

To adapt this script for other OpenFOAM simulations, replace `input_dir` with the
path to your OpenFOAM input files and set the OpenFOAM `distribution` and `version` accordingly.

We run the simulation using the `run` method, specifying the `shell_script` that handles the execution process.

When the simulation is complete, we terminate the machine, download the results and print a summary of the simulation as shown below.

```
Task status: Success

Timeline:
	Waiting for Input         at 21/04, 15:55:50      0.828 s
	In Queue                  at 21/04, 15:55:50      32.128 s
	Preparing to Compute      at 21/04, 15:56:22      3.111 s
	In Progress               at 21/04, 15:56:26      130.295 s
		└> 130.184 s       bash ./Allrun
	Finalizing                at 21/04, 15:58:36      9.277 s
	Success                   at 21/04, 15:58:45      

Data:
	Size of zipped output:    253.64 MB
	Size of unzipped output:  344.54 MB
	Number of output files:   464

Estimated computation cost (US$): 0.0047 US$
```

As you can see in the "In Progress" line, the part of the timeline that represents the actual execution of the simulation, 
the core computation time of this simulation was approximately 2 minutes and 10 seconds.

It's that simple!
