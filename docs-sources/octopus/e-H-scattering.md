# Run the e-H scattering tutorial
This tutorial will show you how to run an official Octopus tutorial using the Inductiva API. 

We will cover the `e-H scattering` use case from the official Octopus tutorials[1], to help you get started with simulations.

## Prerequisites

Start by creating a folder called `SimulationFiles`. Inside that folder you will
need to place the following files:
- [gs.inp](https://storage.googleapis.com/inductiva-api-demo-files/octopus-tutorials/gs.inp)
- [e-h-scat.inp](https://storage.googleapis.com/inductiva-api-demo-files/octopus-tutorials/e-h-scat.inp)
- [plot.gp](https://storage.googleapis.com/inductiva-api-demo-files/octopus-tutorials/plot.gp)

Your `SimulationFiles` should look like this:

```
SimulationFiles/
├── e-h-scat.inp
├── gs.inp
└── plot.gp
```

## Running an Octopus Simulation
Here is the code required to run an Octopus simulation using the Inductiva API:

```python
"""Octopus example"""
import inductiva

# Allocate cloud machine on Google Cloud Platform
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c2d-highcpu-4",
	spot=True)

# Initialize the Simulator
octopus = inductiva.simulators.Octopus( \
    version="16.1")

commands = [
    "mv gs.inp inp",
    "octopus",
    "mv inp gs.inp",
    "mv e-h-scat.inp inp",
    "octopus",
    "mv inp e-h-scat.inp",
    # Generate plots
    "gnuplot plot.gp",
    # Generate gif from plots
    "convert -delay 20 -loop 0 potential_*.png scattering.gif"
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

This example consists of 4 main steps:

1. **Ground-State Calculation:**
   We begin by running `octopus` using the input file `qs.inp`. This step performs the ground-state calculation and generates the `restart/gs` directory, which stores the wavefunctions and other data required for the next stage.

2. **Time-Dependent Simulation:**
   In the second step, we run Octopus with the input file `e-h-scat.inp`. It runs a time-dependent (TD) real-time propagation calculation for a one-dimensional hydrogen-like system using a soft-Coulomb potential. The system starts from scratch and models a single hydrogen atom located at position -10, with an extra electron added, leading to a negatively charged system.
3. **Generating plots**:
   In the third step we run a Gnuplot script that generates time-resolved plots of the electron density and exchange-correlation potential from the Octopus simulation. It creates PNG images at selected time steps, visualizing how the wave packet and potential evolve during the real-time propagation.
4. **Generating gifs**:
   Lastly, we run the command `convert` that will convert all the plots generated into a gif file with an animation of the simulation.

The simulation runs on a cloud machine of type `c2d-highcpu-4`, which provides 4 virtual CPUs. 
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

## Testing different machines

Since this simulation is relatively small, using only 1D and a single atom,
scaling it to a machine with more vCPUs will have limited impact, or none, on
performance. However, to assess potential gains, we will scale from 4 to 8 vCPUs
and test across multiple machine families to compare the performance improvements
offered by different hardware generations.

<div class="cta-bar">
  <div class="cta-text">
    <strong>Start running simulations seamlessly!</strong> You have $5 in <strong>free</strong> credits, no credit card required.
  </div>
  <button  onclick="window.open('https://console.inductiva.ai/', '_blank')" target="_blank" class="cta-button">Sign In</button>
</div>

