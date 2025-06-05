# Run Your First Simulation
This tutorial will show you how to run AMR-Wind simulations using the Inductiva API. 

We will cover the `abl_amd_wenoz` use case from the test files folder of the [AMR-Wind GitHub repository](https://github.com/Exawind/amr-wind/tree/v3.4.0) to help you get started with simulations.

## Prerequisites
Download the required files [here](https://github.com/Exawind/amr-wind/tree/main/test/test_files/abl_amd_wenoz) and place them in a folder called `SimulationFiles`. Then, you’ll be ready to send your simulation to the Cloud.

## Running an AMR-Wind Simulation
Here is the code required to run a AMR-Wind simulation using the Inductiva API:

```python
"""AMR-Wind example."""
import inductiva

# Allocate cloud machine on Google Cloud Platform
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c2d-highcpu-4",
    spot=True)

# Initialize the Simulator
amr_wind = inductiva.simulators.AmrWind(\
    version="3.4.1")

# Run simulation
task = amr_wind.run(input_dir="/Path/to/SimulationFiles",
    sim_config_filename="abl_amd_wenoz.inp",
    on=cloud_machine)

# Wait for the simulation to finish and download the results
task.wait()
cloud_machine.terminate()

task.download_outputs()

task.print_summary()
```

> **Note**: `spot` machines are available at substantial discounts, but your simulation job may be preempted if
> the Cloud provider reclaims the spot machine.

To adapt this script for other AMR-Wind simulations, replace `input_dir` with the
path to your AMR-Wind input files and set the `sim_config_filename` accordingly.

When the simulation is complete, we terminate the machine, download the results and print a summary of the simulation as shown below.

```
Task status: Success

Timeline:
	Waiting for Input         at 16/04, 14:31:00      0.872 s
	In Queue                  at 16/04, 14:31:01      33.07 s
	Preparing to Compute      at 16/04, 14:31:34      1.673 s
	In Progress               at 16/04, 14:31:35      2.954 s
		└> 2.832 s         /opt/openmpi/4.1.6/bin/mpirun --use-hwthread-cpus amr_wind abl_amd_wenoz.inp
	Finalizing                at 16/04, 14:31:38      1.547 s
	Success                   at 16/04, 14:31:40      

Data:
	Size of zipped output:    13.33 MB
	Size of unzipped output:  52.24 MB
	Number of output files:   92

Estimated computation cost (US$): 0.000063 US$
```

As you can see in the "In Progress" line, the part of the timeline that represents the actual execution of the simulation, 
the core computation time of this simulation was approximately 3 seconds.

It's that simple!
