# Run Your First Simulation
This tutorial will show you how to run SPlisHSPlasH simulations using the Inductiva API. 

We will cover the `DamBreakModel` use case from the [SPlisHSPlasH GitHub repository](https://github.com/InteractiveComputerGraphics/SPlisHSPlasH/tree/2.13.0) to help you get started with simulations.

## Prerequisites
1. Download the required `.json` file [here](https://github.com/InteractiveComputerGraphics/SPlisHSPlasH/blob/2.13.0/data/Scenes/DamBreakModel.json) and the additional geometry file [here](https://github.com/InteractiveComputerGraphics/SPlisHSPlasH/blob/2.13.0/data/models/UnitBox.obj). Once downloaded, place both files into a folder named `SimulationFiles`.

2. Make the following adjustments to the `.json` file: update the `geometryFile` path to `UnitBox.obj` and add `"stopAt": 1` as a configuration parameter.

Then, you’ll be ready to send your simulation to the Cloud.

## Running an SPlisHSPlasH Simulation
Here is the code required to run a SPlisHSPlasH simulation using the Inductiva API:

```python
"""SPlisHSPlasH example."""
import inductiva

# Allocate cloud machine on Google Cloud Platform
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c2d-highcpu-4",
    spot=True)

# Initialize the Simulator
splishsplash = inductiva.simulators.SplishSplash(\
    version="2.13.0")

# Run simulation
task = splishsplash.run(input_dir="/Path/to/SimulationFiles",
    sim_config_filename="DamBreakModel.json",
    on=cloud_machine)

# Wait for the simulation to finish and download the results
task.wait()
cloud_machine.terminate()

task.download_outputs()

task.print_summary()
```

> **Note**: `spot` machines are a lot cheaper but may be terminated by the provider if necessary.

To adapt this script for other SPlisHSPlasH simulations, replace `input_dir` with the
path to your SPlisHSPlasH input files and set the `sim_config_filename` accordingly.

When the simulation is complete, we terminate the machine, download the results and print a summary of the simulation as shown below.

```
Task status: Success

Timeline:
	Waiting for Input         at 21/04, 19:49:20      1.116 s
	In Queue                  at 21/04, 19:49:21      34.469 s
	Preparing to Compute      at 21/04, 19:49:55      0.996 s
	In Progress               at 21/04, 19:49:56      9.123 s
		├> 0.982 s         cp /SPlisHSPlasH_CPU/bin/SPHSimulator .
		├> 6.877 s         ./SPHSimulator DamBreakModel.json --no-gui --output-dir .
		└> 1.064 s         rm SPHSimulator
	Finalizing                at 21/04, 19:50:06      0.502 s
	Success                   at 21/04, 19:50:06      

Data:
	Size of zipped output:    1.77 MB
	Size of unzipped output:  6.98 MB
	Number of output files:   5

Estimated computation cost (US$): 0.000099 US$
```

As you can see in the "In Progress" line, the part of the timeline that represents the actual execution of the simulation, 
the core computation time of this simulation was approximately 9 seconds.

It's that simple!
