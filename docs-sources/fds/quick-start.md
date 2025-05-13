# Run Your First Simulation
This tutorial will show you how to run FDS simulations using the Inductiva API. 

We will cover the `box_burn_away6` from the [FDS GitHub repository](https://github.com/firemodels/fds/tree/FDS-6.9.1) to help you get started with simulations.

## Prerequisites
Download the required files [here](https://github.com/firemodels/fds/tree/FDS-6.9.1/Verification/Fires) and place them in a folder called `Fires`. Then, you’ll be ready to send your simulation to the Cloud.

## Running a FDS Simulation
Here is the code required to run a FDS simulation using the Inductiva API:

```python
"""FDS example."""
import inductiva

# Allocate a machine on Google Cloud Platform
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c2d-highcpu-4",
    spot=True)

# Initialize the Simulator
fds = inductiva.simulators.FDS( \
    version="6.9.1")

# Run simulation
task = fds.run( \
    input_dir="/Path/to/Fires",
    sim_config_filename="box_burn_away6.fds",
    n_vcpus=1,
    on=cloud_machine)

# Wait for the simulation to finish and download the results
task.wait()
cloud_machine.terminate()

task.download_outputs()

task.print_summary()
```

> **Note**: `spot` machines are a lot cheaper but may be terminated by the provider if necessary.

Since FDS requires separate mesh setups for each processor, you will need to specify the number of cores (`n_vcpus`) for your simulation. FDS does not automatically assign cores, so it's crucial to configure this manually.

To adapt this script for other FDS simulations, replace `input_dir` with the
path to your FDS input files and set the `sim_config_filename` accordingly.

When the simulation is complete, we terminate the machine, download the results and print a summary of the simulation as shown below.

```
Task status: Success

Timeline:
	Waiting for Input         at 17/04, 15:32:03      0.847 s
	In Queue                  at 17/04, 15:32:04      35.534 s
	Preparing to Compute      at 17/04, 15:32:39      6.411 s
	In Progress               at 17/04, 15:32:46      49.216 s
		└> 49.105 s        /opt/openmpi/4.1.6/bin/mpirun --use-hwthread-cpus --np 1 /opt/fds/Build/ompi_gnu_linux/fds_ompi_gnu_linux box_burn_away6.fds
	Finalizing                at 17/04, 15:33:35      0.697 s
	Success                   at 17/04, 15:33:36      

Data:
	Size of zipped output:    5.64 MB
	Size of unzipped output:  39.30 MB
	Number of output files:   24

Estimated computation cost (US$): 0.00051 US$
```

As you can see in the "In Progress" line, the part of the timeline that represents the actual execution of the simulation, 
the core computation time of this simulation was approximately 49.2 seconds.

It's that simple!
