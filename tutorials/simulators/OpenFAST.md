This guide will walk you through setting up and running OpenFAST simulations using the Inductiva API.

We will cover an example code to help you get started with simulations.

# OpenFAST

[OpenFAST](https://www.nrel.gov/wind/nwtc/openfast.html) provides a 
modular framework for designing and analyzing wind turbines, allowing for 
detailed modeling of components such as blades, towers, and control systems.

OpenFAST comes integrated with FAST.Farm, an extension designed for 
large-scale wind farm simulations. FAST.Farm models complex interactions 
between turbines, accounting for factors like wake effects, turbulence, 
and terrain. FAST.Farm's capabilities are indispensable for the design 
and planning of wind energy projects, unlocking the full potential of 
wind resources for renewable energy generation.

Both OpenFAST and FAST.Farm are compiled with OpenMP to enable parallel processing, enhancing the performance and scalability of your simulations.

### Objective

We will cover a quick start use case available in the [OpenFAST GitHub Repository](https://github.com/openfast) to help you get started with simulations.

This use case uses the NREL 5-MW wind turbine, a hypothetical yet representative multi-MW wind turbine with a rated power of 5 MW, a rated rotor speed of 12.1 rpm, a hub height of 90 m, and a rotor diameter of 126 m. It focuses on an “onshore” version
of the turbine, solely considering the structure (no aerodynamics), where the tower is initially offset by 3 m at the top.

### Prerequisites

Download the required files [here](https://github.com/OpenFAST/r-test/tree/main/glue-codes/openfast/MinimalExample) and place them in a folder called `MinimalExample`. Then, you’ll be ready to send your simulation to the Cloud.

### Running Your Simulation

Here is the code required to run an OpenFAST simulation using the Inductiva API:

```python
"""OpenFAST example."""
import inductiva

# Allocate cloud machine on Google Cloud Platform
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c2-standard-4")

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

When the simulation is complete, we terminate the machine, download the results and print a summary of the simulation as shown below.

```
Task status: Success

Timeline:
	Waiting for Input         at 10/03, 20:20:43      0.935 s
	In Queue                  at 10/03, 20:20:44      30.778 s
	Preparing to Compute      at 10/03, 20:21:15      1.353 s
	In Progress               at 10/03, 20:21:16      3.89 s
		└> 3.773 s         openfast Main.fst
	Finalizing                at 10/03, 20:21:20      0.408 s
	Success                   at 10/03, 20:21:20      

Data:
	Size of zipped output:    47.76 KB
	Size of unzipped output:  141.94 KB
	Number of output files:   4

Estimated computation cost (US$): 0.00011 US$
```

As you can see in the "In Progress" line, the part of the timeline that represents the actual execution of the simulation, 
the core computation time of this simulation was 3.8 seconds.

It's that simple!