# Run Your First Simulation
 A typical GROMACS simulation involves several key steps: preparing the molecular structure, performing energy minimization, running the simulation, and analyzing the results. You’ll need multiple input files to configure and run your simulation, and specific GROMACS commands will be used to carry out each step.

This tutorial will show you how to run GROMACS simulations using the Inductiva API. 

We will cover the `Bulk salt solution` use case from the [GROMACS tutorials page](https://gromacstutorials.github.io/sphinx/build/html/tutorials/tutorial1/bulk-solution.html) to help you get started with simulations.

## Prerequisites
To generate the input files required for this simulation, please complete the initial steps of [this section](https://gromacstutorials.github.io/sphinx/build/html/tutorials/tutorial1/bulk-solution.html) of the tutorial, as well as the `Energy Minimization` part.

After completing these steps, you should have the following files in your `SimulationFiles` folder:

```
├── ./conf.gro
├── ./ff
│   ├── ./ff/forcefield.itp
│   ├── ./ff/h2o.itp
│   ├── ./ff/na.itp
│   └── ./ff/so4.itp
├── ./min.mdp
└── ./topol.top
```

With these files in place, you're ready to send your simulation to the cloud!

## Running an GROMACS Simulation
Here is the code required to run a GROMACS simulation using the Inductiva API:

```python
"""GROMACS example."""
import inductiva

# Allocate cloud machine on Google Cloud Platform
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c2d-highcpu-16",
    spot=True)

# List of commands to run
commands = [
    "gmx grompp -f min.mdp -c conf.gro -p topol.top -o min -pp min -po min",
    "gmx mdrun -v -deffnm min",
]

# Initialize the Simulator
gromacs = inductiva.simulators.GROMACS( \
    version="2025.0")

# Run simulation
task = gromacs.run( \
    input_dir="/Path/to/SimulationFiles",
    commands=commands,
    on=cloud_machine)

# Wait for the simulation to finish and download the results
task.wait()
cloud_machine.terminate()

task.download_outputs()

task.print_summary()
```

> **Note**: `spot` machines are a lot cheaper but may be terminated by the provider if necessary.

To adapt this script for other GROMACS simulations, replace `input_dir` with the
path to your GROMACS input files and set the `commands` accordingly.

When the simulation is complete, we terminate the machine, download the results
and print a summary of the simulation as shown below.

```
Task status: Success

Timeline:
    Waiting for Input         at 22/04, 09:28:14      0.875 s
    In Queue                  at 22/04, 09:28:15      34.39 s
    Preparing to Compute      at 22/04, 09:28:49      1.699 s
    In Progress               at 22/04, 09:28:51      3.056 s
        ├> 1.915 s         gmx grompp -f min.mdp -c conf.gro -p topol.top -o min -pp min -po min
        └> 0.983 s         gmx mdrun -v -deffnm min
    Finalizing                at 22/04, 09:28:54      0.485 s
    Success                   at 22/04, 09:28:54      

Data:
    Size of zipped output:    1.75 MB
    Size of unzipped output:  2.29 MB
    Number of output files:   9

Estimated computation cost (US$): 0.00020 US$
```

As you can see in the "In Progress" line, the part of the timeline that
represents the actual execution of the simulation, the core computation time of
this simulation was approximately 3.1 seconds.

It's that simple!
