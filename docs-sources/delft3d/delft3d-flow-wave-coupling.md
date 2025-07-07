# Run Delft3D with FLOW-WAVE coupling
This tutorial demonstrates how to run a Delft3D simulation with FLOW-WAVE coupling using the Inductiva API.

We’ll guide you through the `03_flow-wave example` provided in the [Delft3D Subversion repository](https://svn.oss.deltares.nl/repos/delft3d/branches/releases/7545/).

## Prerequisites
Before running the simulation, you'll need to download the input files for the `03_flow-wave` example.

### Option 1: Download manually
You can download the files directly from the [Delft3D repository](https://svn.oss.deltares.nl/repos/delft3d/branches/releases/7545/examples/03_flow-wave/) and save them in a folder named `03_flow-wave`.

### Option 2: Use SVN
If you have Subversion (SVN) installed, you can quickly download all files by running this command in your terminal:

```
svn checkout https://svn.oss.deltares.nl/repos/delft3d/branches/releases/7545/examples/03_flow-wave/
```

This will create a folder called `03_flow-wave` containing all the required files.

## Creating the Simulation Script
This type of simulation is slightly more complex, as it involves running two components concurrently:
- `d_hydro.exe` for the flow module
- `wave.exe` for the wave module

To manage this, you’ll create a shell script that launches both commands simultaneously.

In the `03_flow-wave` directory, create a file named `run_sim.sh` with the following content:

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

The `D3D_HOME` and `ARCH` environment variables are automatically set by the Delft3D environment, so no manual configuration is needed.

Once you have your executables defined, you just need to run the `d_hydro.exe` and `wave.exe` commands with the appropriate arguments. In this case, we use `mpirun` to run `d_hydro.exe` with the specified number of processes (`procs`), followed by `wave.exe` using the specified input file (`mdwfile`).

> **Note**: The `d_hydro.exe` command ends with an `&`, which runs it in the background. This allows `wave.exe` to run in the foreground. Running both simultaneously is necessary, as they are coupled and need to operate together throughout the simulation.

## Running the Coupled Simulation
Here is the code required to run this Delft3D simulation using the Inductiva API:

```python
"""Delft3D FLOW-WAVE coupled simulation."""
import inductiva

# Allocate a machine on Google Cloud Platform
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c2d-highcpu-16",
	spot=True)

# Initialize the Simulator
delft3d = inductiva.simulators.Delft3D(\
    version="6.04.00")

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
```

As you can see in the "In Progress" line, the part of the timeline that represents the actual execution of 
the simulation, the core computation time of this simulation was approximately 85.2 seconds.

And that’s it — you’re now ready to run coupled FLOW-WAVE simulations using Delft3D on the cloud with Inductiva!
