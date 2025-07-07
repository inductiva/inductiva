# Run Delft3D with FLOW-WAVE coupling
This tutorial will show you how to run Deft3D with FLOW-WAVE coupling using the Inductiva API. 

We will cover the `03_flow-wave` use case from the examples available in the [Delft3D Subversion repository](https://svn.oss.deltares.nl/repos/delft3d/branches/releases/7545/).

## Prerequisites
Download the required files [here](https://svn.oss.deltares.nl/repos/delft3d/branches/releases/7545/examples/03_flow-wave/) and save them to a folder named `03_flow-wave`.

**Notes:** To download the files run `svn checkout https://svn.oss.deltares.nl/repos/delft3d/branches/releases/7545/examples/03_flow-wave/` on your terminal.

## Creating Your Simulation Script

This type of simulation is a bit more complex than the previous ones, as it
requires running two components simultaneously: `d_hydro` for the flow module
and `wave.exe` for the wave module. To manage this, we’ll create a shell script
that launches both commands together.

Create a file named `run_sim.sh` in the `03_flow-wave` directory with the following content:

```bash
#!/bin/bash

argfile=config_d_hydro.xml
mdwfile=r17.mdw
procs=16

flowexedir=$D3D_HOME/$ARCH/flow2d3d/bin
waveexedir=$D3D_HOME/$ARCH/wave/bin
swanexedir=$D3D_HOME/$ARCH/swan/bin
swanbatdir=$D3D_HOME/$ARCH/swan/scripts

# Run
# Start FLOW
mpirun -np $procs $flowexedir/d_hydro.exe $argfile &

# Start WAVE
$waveexedir/wave.exe $mdwfile 1
```

Notice that we’re referencing the required executables using the `D3D_HOME` and
`ARCH` environment variables. These are set automatically by the Delft3D
environment, so there’s no need to define them manually.

Once you have your executables defined you just need to run the `d_hydro.exe` and `wave.exe` commands
with the appropriate arguments. In this case, we are using `mpirun` to run
`d_hydro.exe` with the specified number of processes (`procs`), and then running
`wave.exe` with the specified input file (`mdwfile`).

**Notes:** See that we ran `d_hydro.exe` with `&` at the end of the command. This is to run it in the background, so that we can run `wave.exe` in the foreground. This is necessary because `d_hydro.exe` and `wave.exe` will run both coupled with each other.

## Running the coupled simulation
Here is the code required to run the Delft3D simulation using the Inductiva API:

```python
"""Delft3D Simulation."""
import inductiva

# Allocate a machine on Google Cloud Platform
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c2d-highcpu-16",
	spot=True)

# Initialize the Simulator
delft3d = inductiva.simulators.Delft3D(\
    version="6.04.00",use_dev=True)

# Run simulation
task = delft3d.run( \
    input_dir="/Path/to/03_flow-wave",
    commands = ["bash run_sim.sh"],
    on=cloud_machine)

# Wait for the simulation to finish and download the results
task.wait()
cloud_machine.terminate()

task.download_outputs()

task.print_summary()
```

In this basic example, we're using a cloud machine (`c2d-highcpu-16`) equipped with 16 virtual CPUs. 
For larger or more compute-intensive simulations, consider adjusting the `machine_type` parameter to select 
a machine with more virtual CPUs and increased memory capacity. You can explore the full range of available machines [here](https://console.inductiva.ai/machine-groups/instance-types).

> **Note**: Setting `spot=True` enables the use of spot machines, which are available at substantial discounts. 
> However, your simulation may be interrupted if the cloud provider reclaims the machine.

To adapt this script for other Delft3D simulations, replace `input_dir` with the
path to your Delft3D input files and set the the `commands` accordingly.

When the simulation is complete, we terminate the machine, download the results and print a summary of the simulation as shown below.

```
Task status: Success

Timeline:
	Waiting for Input         at 03/07, 11:45:29      1.101 s
	In Queue                  at 03/07, 11:45:30      46.132 s
	Preparing to Compute      at 03/07, 11:46:16      2.973 s
	In Progress               at 03/07, 11:46:19      85.204 s
		└> 85.022 s        bash run_sim.sh
	Finalizing                at 03/07, 11:47:44      0.504 s
	Success                   at 03/07, 11:47:45      

Data:
	Size of zipped output:    5.52 MB
	Size of unzipped output:  18.39 MB
	Number of output files:   76

Estimated computation cost (US$): 0.0023 US$

Go to https://console.inductiva.ai/tasks/975yio4uujcz7hssjszbdupdi for more details.
```

As you can see in the "In Progress" line, the part of the timeline that represents the actual execution of the simulation, 
the core computation time of this simulation was approximately 85.2 seconds.

It's that simple!
