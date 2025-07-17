# Run the e-H scattering tutorial
This tutorial will show you how to run an official Octopus tutorial using the Inductiva API. 

We will cover the `e-H scattering` use case from the official Octopus tutorials[1], to help you get started with simulations.

<p align="center">
  <img src="./_static/scattering.gif" alt="e-H Simulation" width="50%" height="50%">
</p>


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
	Waiting for Input         at 15/07, 11:21:09      0.743 s
	In Queue                  at 15/07, 11:21:10      41.361 s
	Preparing to Compute      at 15/07, 11:21:52      4.081 s
	In Progress               at 15/07, 11:21:56      685.099 s
		├> 1.148 s         mv gs.inp inp
		├> 2.087 s         octopus
		├> 1.084 s         mv inp gs.inp
		├> 1.092 s         mv e-h-scat.inp inp
		├> 472.522 s       octopus
		├> 1.074 s         mv inp e-h-scat.inp
		├> 7.09 s          gnuplot plot.gp
		└> 198.277 s       convert -delay 20 -loop 0 potential_*.png scattering.gif
	Finalizing                at 15/07, 11:33:21      1.349 s
	Success                   at 15/07, 11:33:22      

Data:
	Size of zipped output:    38.64 MB
	Size of unzipped output:  102.43 MB
	Number of output files:   2777

Estimated computation cost (US$): 0.0045 US$
```

As you can see in the "In Progress" line, the part of the timeline that
represents the actual execution of the simulation, the core computation time of
this simulation was approximately 11 minutes and 25 seconds.

## Testing Different Machines

This simulation is relatively small, using only a one-dimensional setup and a
single atom. As a result, increasing the number of vCPUs is expected to have
minimal or no impact on performance. Nevertheless, to explore potential speedups,
we tested configurations with 4 and 8 vCPUs across several machine families to
compare performance across hardware generations.

| Machine Type  | Execution Time         | Estimated Cost |
| ------------- | ---------------------- | -------------- |
| c2d-highcpu-4 | 11 minutes, 25 seconds | 0.0045 US$    |
| c2d-highcpu-8 | 11 minutes, 22 seconds | 0.0082 US$    |
| c3d-highcpu-4 | 11 minutes, 15 seconds | 0.0079 US$    |
| c3d-highcpu-8 | 11 minutes, 11 seconds | 0.0150 US$    |
| c4-highcpu-4  | 6 minutes, 51 seconds  | 0.0092 US$    |
| c4-highcpu-8  | 7 minutes, 14 seconds  | 0.0190 US$    |

The results confirm that increasing the number of vCPUs has little effect on
performance for such a small simulation, with runtimes remaining nearly identical
between 4 and 8 vCPUs. The newer `c4-highcpu` machines complete
the task much faster than the `c2d` and `c3d` families, though at a higher cost.
Among all options, `c2d-highcpu-4` is the most cost-effective, while `c4-highcpu-4`
offers the best performance. This highlights the limited benefit of scaling CPU
count for this particular simulation.

<div class="cta-bar">
  <div class="cta-text">
    <strong>Start running simulations seamlessly!</strong> You have $5 in <strong>free</strong> credits, no credit card required.
  </div>
  <button  onclick="window.open('https://console.inductiva.ai/', '_blank')" target="_blank" class="cta-button">Sign In</button>
</div>

**References:**

[[1]Official Octopus tutorial](https://octopus-code.org/documentation/13/tutorial/model/e-h_scattering/)