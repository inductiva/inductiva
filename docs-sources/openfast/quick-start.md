# Run Your First Simulation
This tutorial will show you how to run OpenFAST simulations using the Inductiva API. 

We will cover a quick start use case available in the [OpenFAST GitHub repository](https://github.com/openfast) to help you get started with simulations.

This use case uses the NREL 5-MW wind turbine, a hypothetical yet representative multi-MW wind turbine with a rated power of 5 MW, a rated rotor speed of 12.1 rpm, 
a hub height of 90 m, and a rotor diameter of 126 m. It focuses on an “onshore” version of the turbine, considering only the structure (no aerodynamics), where the top of the tower is initially displaced horizontally by 3 m from its equilibrium position.

## Prerequisites
Download the required files [here](https://github.com/OpenFAST/r-test/tree/main/glue-codes/openfast/MinimalExample) and place them in a folder called `MinimalExample`. Then, you’ll be ready to send your simulation to the Cloud.

## Running an OpenFAST Simulation
Here is the code required to run an OpenFAST simulation using the Inductiva API:

```python
"""OpenFAST example."""
import inductiva

# Allocate cloud machine on Google Cloud Platform
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c2d-highcpu-4",
	spot=True)

# List of commands to run
commands = ["openfast Main.fst"]

# Initialize the Simulator
openfast = inductiva.simulators.OpenFAST( \
    version="4.0.2")

# Run simulation
task = openfast.run(input_dir="/Path/to/MinimalExample",
              sim_config_filename="Main.fst",
              commands=commands,
              on=cloud_machine)

# Wait for the simulation to finish and download the results
task.wait()
cloud_machine.terminate()

task.download_outputs()

task.print_summary()
```

> **Note**: `spot` machines are a lot cheaper but may be terminated by the provider if necessary.

To adapt this script for other OpenFAST simulations, replace `input_dir` with the
path to your OpenFAST input files and set the `sim_config_filename` accordingly.

When the simulation is complete, we terminate the machine, download the results and print a summary 
of the simulation as shown below.

```
Task status: Success

Timeline:
	Waiting for Input         at 21/04, 15:44:16      0.959 s
	In Queue                  at 21/04, 15:44:17      31.921 s
	Preparing to Compute      at 21/04, 15:44:49      1.346 s
	In Progress               at 21/04, 15:44:50      3.45 s
		└> 3.336 s         openfast Main.fst
	Finalizing                at 21/04, 15:44:53      0.4 s
	Success                   at 21/04, 15:44:54      

Data:
	Size of zipped output:    47.58 KB
	Size of unzipped output:  141.94 KB
	Number of output files:   3

Estimated computation cost (US$): 0.000054 US$
```

As you can see in the "In Progress" line, the part of the timeline that represents the actual execution of the simulation, 
the core computation time of this simulation was approximately 3.5 seconds.

It's that simple!
