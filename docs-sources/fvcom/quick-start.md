# Run Your First Simulation
This tutorial will show you how to run FVCOM simulations using the Inductiva API. 

We will cover the `Estuary` use case from the test files folder of the [FVCOM GitHub repository](https://github.com/FVCOM-GitHub/FVCOM) to help you get started with simulations.

## Prerequisites
Download the required files [here](https://storage.googleapis.com/inductiva-api-demo-files/fvcom-input-example.zip). Then, you’ll be ready to send your simulation to the Cloud.

## Running a FVCOM Simulation
Here is the code required to run a FVCOM simulation using the Inductiva API:

```python
"""FVCOM Simulation."""
import inductiva

# Allocate cloud machine on Google Cloud Platform
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c2d-highcpu-4",
	spot=True)

# Initialize the Simulator
fvcom = inductiva.simulators.FVCOM(\
    version="5.1.0")

# Run simulation
task = fvcom.run( \
    input_dir="/Path/to/fvcom-input-example",
    working_dir="run/",
    case_name="tst",
    on=cloud_machine)

# Wait for the simulation to finish and download the results
task.wait()
cloud_machine.terminate()

task.download_outputs()

task.print_summary()
```

> **Note**: `spot` machines are a lot cheaper but may be terminated by the provider if necessary.

To adapt this script for other FVCOM simulations, replace `input_dir` with the
path to your FVCOM input files and set the `case_name` accordingly.

When the simulation is complete, we terminate the machine, download the results and print a summary of the simulation as shown below.

```
Task status: Success

Timeline:
	Waiting for Input         at 17/04, 15:39:20      0.879 s
	In Queue                  at 17/04, 15:39:21      38.245 s
	Preparing to Compute      at 17/04, 15:39:59      1.686 s
	In Progress               at 17/04, 15:40:01      2.029 s
		└> 1.904 s         /opt/openmpi/4.1.6/bin/mpirun --use-hwthread-cpus fvcom --CASENAME=tst --dbg=0
	Finalizing                at 17/04, 15:40:03      0.42 s
	Success                   at 17/04, 15:40:03      

Data:
	Size of zipped output:    72.36 KB
	Size of unzipped output:  1.72 MB
	Number of output files:   5

Estimated computation cost (US$): 0.000045 US$
```

As you can see in the "In Progress" line, the part of the timeline that represents the actual execution of the simulation, 
the core computation time of this simulation was approximately 2 seconds.

It's that simple!
