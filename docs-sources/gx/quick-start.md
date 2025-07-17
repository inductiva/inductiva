# Run Your First Simulation
This tutorial will show you how to run GX simulations using the Inductiva API. 

This tutorial will cover a non-linear use case, available in the [GX documentation](https://gx.readthedocs.io/en/latest/Nonlinear.html), to help you get started with simulations.

## Prerequisites
Download the required files [here](https://bitbucket.org/gyrokinetics/gx/src/gx/benchmarks/nonlinear/cyclone/cyclone_miller_adiabatic_electrons.in) and place them in a folder called `NonlinearExample`. Then, you’ll be ready to send your simulation to the Cloud.

## Running a GX Simulation
Here is the code required to run GX simulation using the Inductiva API:

```python
"""GX Simulation"""
import inductiva

# Allocate a machine on Google Cloud Platform
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="g2-standard-4",
	spot=True)

# Initialize the Simulator
gx = inductiva.simulators.GX(\
    version="v.11-2024")

# Run simulation
task = gx.run( \
    input_dir="/Path/to/NonlinearExample",
    sim_config_filename="cyclone_miller_adiabatic_electrons.in",
    on=cloud_machine)

# Wait for the simulation to finish and download the results
task.wait()
cloud_machine.terminate()

task.download_outputs()

task.print_summary()
```

In this basic example, we're using a cloud machine (`g2-standard-4`) equipped with 4 virtual CPUs and an NVIDIA L4 GPU. 
For larger or more compute-intensive simulations, consider adjusting the `machine_type` parameter to select 
a more powerful GPU machine. You can explore the full range of available 
machines [here](https://console.inductiva.ai/machine-groups/instance-types).

> **Note**: Setting `spot=True` enables the use of spot machines, which are available at substantial discounts. 
> However, your simulation may be interrupted if the cloud provider reclaims the machine.

To adapt this script for other GX simulations, replace `input_dir` with the
path to your GX input files and set the `sim_config_filename` accordingly.

When the simulation is complete, we terminate the machine, download the results and print a summary of the simulation as shown below.

```
Task status: Success

Timeline:
	Waiting for Input         at 28/03, 14:26:20      0.787 s
	In Queue                  at 28/03, 14:26:21      31.739 s
	Preparing to Compute      at 28/03, 14:26:53      20.526 s
	In Progress               at 28/03, 14:27:14      4.36 s
		└> 4.064 s         gx cyclone_miller_adiabatic_electrons.in
	Finalizing                at 28/03, 14:27:18      0.432 s
	Success                   at 28/03, 14:27:18      

Data:
	Size of zipped output:    1.63 KB
	Size of unzipped output:  9.05 KB
	Number of output files:   4

Estimated computation cost (US$): 0.0022 US$
```

As you can see in the "In Progress" line, the part of the timeline that represents the actual execution of the simulation, 
the core computation time of this simulation was approximately 4 seconds.

<div class="cta-bar">
  <div class="cta-text">
  	<strong>Kickstart your simulations!</strong> You have $5 in <strong>free credits</strong>, no credit card required.
  </div>
 <button  onclick="window.open('https://console.inductiva.ai/?utm_source=guide_gx&utm_medium=button&utm_campaign=signup', '_blank')" target="_blank" class="cta-button">Sign In</button>
</div>

