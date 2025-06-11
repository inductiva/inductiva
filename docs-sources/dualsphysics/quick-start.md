# Run Your First Simulation

## Technical Details
Running your DualSPHysics simulation workflows on Inductiva is similar to running them locally, 
with a few key differences. Instead of executing your DualSPHysics shell script directly, 
you'll need to pass it to the Inductiva API using the run() method. This will allow your simulation to be 
executed on a remote resource. If you already have a functioning shell script that orchestrates your DualSPHysics 
simulation locally, you're almost ready to run it on Inductiva.

### Key Adaptations for Inductiva
While the overall process is similar, there are a few minor adjustments that might be needed to accommodate 
the different environment:

- **Path-related changes**: The location of the DualSPHysics binaries in the Inductiva infrastructure may differ from your local setup. You may need to modify path variables in your orchestration script.
- **Interactive commands**: Any commands that require keyboard input will need to be adjusted. You can either remove these commands or set them to run with default parameters.

### Location of binaries
All DualSPHysics binaries are located in the directory `/DualSPHysics_v5.2/bin/linux/`.
This path is included in the system's `PATH` environment variable, so you
can execute commands directly without needing to specify the full path.

We've also created symbolic links for easier access. For example, instead
of using the full name `DualSPHysics5.2CPU_linux64`, you can simply call
`dualsphysics`. The naming convention is simple: we remove the `_linux64`
suffix and convert the name to lowercase.

This renaming allows you to use either the full binary name (e.g.,
`DualSPHysics5.2CPU_linux64`) or the simplified name (e.g., `dualsphysics`),
which abstracts away the architecture and simulator version.

## Objective
This tutorial will show you how to run DualSPHysics simulations using the Inductiva API. 

We will cover a classical CFD case of a flow over a cylinder from the [DualSPHysics offical download page](https://dual.sphysics.org/downloads/) to help you get started with simulations.

## Prerequisites
Download the required files present in `DualSPHysics_v5.2.2.zip`
[here](https://dual.sphysics.org/downloads/) and place them in a folder called
`DualSPHysics_v5.2`. All your simulation files will be located at `DualSPHysics_v5.2/examples/inletoutlet/01_FlowCylinder`.

After that you will need to make some changes to the file `xCaseFlowCylinder_Re020_linux64_GPU.sh`:

- Update the `dirbin` to `export dirbin=/DualSPHysics_v5.2/bin/linux/`
- Remove the last line of the script `read -n1 -r -p "Press any key to continue..." key`

You now have all you need to start your simulation at `DualSPHysics_v5.2/examples/inletoutlet/01_FlowCylinder`.

## Running a DualSPHysics Simulation
Here is the code required to run a DualSPHysics simulation using the Inductiva API:

```python
"""DualSPHysics example."""
import inductiva

# Allocate a machine on Google Cloud Platform
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="g2-standard-32",
    spot=True)

# Initialize the Simulator
dualsphysics = inductiva.simulators.DualSPHysics( \
    version="5.2.1")

# Run simulation
task = dualsphysics.run( \
    input_dir="/Path/to/01_FlowCylinder",
    shell_script="xCaseFlowCylinder_Re020_linux64_GPU.sh",
    on=cloud_machine)

# Wait for the simulation to finish and download the results
task.wait()
cloud_machine.terminate()

task.download_outputs()

task.print_summary()
```

In this basic example, we're using a cloud machine (`g2-standard-32`) equipped with 32 virtual CPUs and an NVIDIA L4 GPU. 
For larger or more compute-intensive simulations, consider adjusting the `machine_type` parameter to select 
a more powerful GPU machine. You can explore the full range of available 
machines [here](https://console.inductiva.ai/machine-groups/instance-types).

> **Note**: Setting `spot=True` enables the use of spot machines, which are available at substantial discounts. 
> However, your simulation may be interrupted if the cloud provider reclaims the machine.

To adapt this script for other DualSPHysics simulations, replace `input_dir` with the
path to your DualSPHysics input files and set the `shell_script` accordingly.

When the simulation is complete, we terminate the machine, download the results and print a summary of the simulation as shown below.

```
Task status: Success

Timeline:
	Waiting for Input         at 09/04, 10:29:35      0.871 s
	In Queue                  at 09/04, 10:29:36      54.536 s
	Preparing to Compute      at 09/04, 10:30:31      12.672 s
	In Progress               at 09/04, 10:30:44      99.348 s
		└> 99.186 s        bash xCaseFlowCylinder_Re020_linux64_GPU.sh
	Finalizing                at 09/04, 10:32:23      64.657 s
	Success                   at 09/04, 10:33:28      

Data:
	Size of zipped output:    1.25 GB
	Size of unzipped output:  2.41 GB
	Number of output files:   2032

Estimated computation cost (US$): 0.034 US$

Go to https://console.inductiva.ai/tasks/pn5bnuj7ygsenouwnm7zmxr2s for more details.
```

As you can see in the "In Progress" line, the part of the timeline that
represents the actual execution of the simulation, 
the core computation time of this simulation was 99.3 seconds (around 1 minute
and 39 seconds).

It's that simple!
