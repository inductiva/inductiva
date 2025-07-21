# Run Your First Simulation
This tutorial will show you how to run Octopus simulations using the Inductiva API. 

We will cover the `01-propagators.03-etrs_taylor` use case from the Octopus Gitlab repository[1], to help you get started with simulations.

## Prerequisites
Download the following files and save them into a folder called `SimulationFiles` with the
respective name:
- Download [this](https://gitlab.com/octopus-code/octopus/-/raw/16.1/testsuite/real_time/01-propagators.03-etrs_taylor.inp?ref_type=tags) file and save it as `td.inp`
- Download [this](https://gitlab.com/octopus-code/octopus/-/raw/16.1/testsuite/real_time/02-propagators.01-gs.inp?ref_type=tags) file and save it as `gs.inp`

Your `SimulationFiles` folder should look like this:

```
SimulationFiles/
├── gs.inp
└── td.inp
```

Now, you are ready to send your simulation to the Cloud.

## Running an Octopus Simulation
Here is the code required to run an Octopus simulation using the Inductiva API:

```python
"""Octopus example"""
import inductiva

# Allocate cloud machine on Google Cloud Platform
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c2d-highcpu-16",
	spot=True)

# Initialize the Simulator
octopus = inductiva.simulators.Octopus( \
    version="16.1",
    use_dev=True)

commands = [
    "mv gs.inp inp",
    "octopus",
    "mv inp gs.inp",
    "mv td.inp inp",
    "octopus",
    "mv inp td.inp"
]

# Run simulation
task = octopus.run( \
    input_dir="/Path/to/SimulationFiles",
    commands=commands,
    on=cloud_machine)

task.wait()
cloud_machine.terminate()

task.download_outputs()

task.print_summary()
```

This example consists of two main steps:

1. **Ground-State Calculation:**
   We begin by running `octopus` using the input file `qs.inp`. This step performs the ground-state calculation and generates the `restart/gs` directory, which stores the wavefunctions and other data required for the next stage.

2. **Time-Dependent Simulation:**
   In the second step, we run Octopus with the input file `td.inp`. This launches a real-time, time-dependent simulation based on **Time-Dependent Density Functional Theory (TDDFT)**, using the results from the previous ground-state step.

The simulation runs on a cloud machine of type `c2d-highcpu-16`, which provides 16 virtual CPUs. 
For larger or more compute-intensive simulations, consider adjusting the `machine_type` parameter to select 
a machine with more virtual CPUs and increased memory capacity. You can explore the full range of available machines [here](https://console.inductiva.ai/machine-groups/instance-types).

> **Note**: Setting `spot=True` enables the use of spot machines, which are available at substantial discounts. 
> However, your simulation may be interrupted if the cloud provider reclaims the machine.

To adapt this script for other Octopus simulations, replace `input_dir` with the
path to your Octopus input files and set the Octopus `version` accordingly.

When the simulation is complete, we terminate the machine, download the results and print a summary of the simulation as shown below.

```
Task status: Success

Timeline:
	Waiting for Input         at 14/07, 15:29:45      0.771 s
	In Queue                  at 14/07, 15:29:46      49.68 s
	Preparing to Compute      at 14/07, 15:30:36      5.735 s
	In Progress               at 14/07, 15:30:42      9.018 s
		├> 1.086 s         mv gs.inp inp
		├> 2.093 s         octopus
		├> 1.083 s         mv inp gs.inp
		├> 1.073 s         mv td.inp inp
		├> 2.083 s         octopus
		└> 1.075 s         mv inp td.inp
	Finalizing                at 14/07, 15:30:51      0.529 s
	Success                   at 14/07, 15:30:51      

Data:
	Size of zipped output:    557.60 KB
	Size of unzipped output:  1.49 MB
	Number of output files:   51

Estimated computation cost (US$): 0.00039 US$
```

As you can see in the "In Progress" line, the part of the timeline that represents the actual execution of the simulation, 
the core computation time of this simulation was approximately 9 seconds.

```{banner_small}
:origin: octopus
```

**References:**

[[1]Octopus Gitlab repository](https://gitlab.com/octopus-code/octopus/-/blob/16.1/testsuite/real_time/01-propagators.03-etrs_taylor.inp?ref_type=tags)