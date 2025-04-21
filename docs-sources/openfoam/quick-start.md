# Run Your First Simulation
This tutorial will show you how to run OpenFOAM simulations using the Inductiva API. 

We will cover the `motorbike` use case from the the [OpenFOAM Foundation GitHub repository](https://github.com/OpenFOAM/OpenFOAM-8/tree/version-8/tutorials), to help you get started with simulations.

## Prerequisites
Download the required files [here](https://github.com/OpenFOAM/OpenFOAM-8/tree/version-8/tutorials/incompressible/simpleFoam/motorBike) and place them in a folder called `SimulationFiles`. Then, you’ll be ready to send your simulation to the Cloud.

## Running an OpenFOAM Simulation
Here is the code required to run a OpenFOAM simulation using the Inductiva API:

```python
"""OpenFOAM example"""
import inductiva

# Allocate cloud machine on Google Cloud Platform
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c3d-standard-16",
	spot=True)

# Initialize the Simulator
OpenFOAM = inductiva.simulators.OpenFOAM( \
    version="8",
	distribution="foundation")

# Run simulation
task = OpenFOAM.run(input_dir="/Path/to/RegularWavePropagation",
    shell_script="./Allrun",
    on=cloud_machine)

# Wait for the simulation to finish and download the results
task.wait()
cloud_machine.terminate()

task.download_outputs()

task.print_summary()
```

> **Note**: `spot` machines are a lot cheaper but may be terminated by the provider if necessary.

To adapt this script for other OpenFOAM simulations, replace `input_dir` with the
path to your OpenFOAM input files and set the OpenFOAM `distribution` and `version` accordingly.

We run the simulation using the `run` method, specifying the `shell_script` that handles the execution process.

When the simulation is complete, we terminate the machine, download the results and print a summary of the simulation as shown below.

```
Task status: Success

Timeline:
	Waiting for Input         at 10/04, 10:44:12      0.78 s
	In Queue                  at 10/04, 10:44:12      40.014 s
	Preparing to Compute      at 10/04, 10:44:52      3.416 s
	In Progress               at 10/04, 10:44:56      185.345 s
		└> 185.228 s       bash ./Allrun
	Finalizing                at 10/04, 10:48:01      14.257 s
	Success                   at 10/04, 10:48:15      

Data:
	Size of zipped output:    253.64 MB
	Size of unzipped output:  344.54 MB
	Number of output files:   463

Estimated computation cost (US$): 0.011 US$
```

As you can see in the "In Progress" line, the part of the timeline that represents the actual execution of the simulation, 
the core computation time of this simulation was approximately 3 minutes.

It's that simple!
