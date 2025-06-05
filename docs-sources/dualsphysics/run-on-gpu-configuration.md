# Run DualSPHysics Across Different GPUs
In this tutorial, we’ll demonstrate how to use the Inductiva API to run DualSPHysics simulations on various GPU configurations.

To illustrate this, we’ll run a benchmark simulation based on the study presented in the paper
[A fluid–structure interaction model for free-surface flows and flexible structures using smoothed particle hydrodynamics on a GPU](https://www.sciencedirect.com/science/article/pii/S0889974621000955?via%3Dihub).

Among the simulations discussed, we focus on the **3-D dam break impacting an elastic plate** scenario.

<p align="center"><img src="./_static/dam_break_elastic.gif" alt="Visualization created with Blender." width="700"></p>

## Simulation Overview
This simulation illustrates the classic dam break problem, where a column of water is suddenly released. What makes it unique is the presence of an elastic plate positioned at the opposite side of the domain, interacting dynamically with the fluid flow.

> For more details, see section 4.5, *3-D dam break impacting an elastic plate*, in the referenced paper.

## Prerequisites

### Download the DualSPHysics package
Download the required files from `DualSPHysics_v5.4.2.zip` [here](https://dual.sphysics.org/downloads/).
You will be working within this directory and writing the Inductiva Python script there.

### Update the simulation script of the `examples/flexstruc/01_DamBreak` case

Before running the simulation, update the script located at: `examples/flexstruc/01_DamBreak/xCaseDambreak3D_FSI_linux64_GPU.sh`.

Make the following adjustments:
1. **Update the `dirbin` variable:**
   Modify the `xCaseDambreak3D_FSI_linux64_GPU.sh` script to point to the correct binaries directory:
   ```bash
   export dirbin=/DualSPHysics/bin/linux/
   ```
2. **Remove user input prompt:**
   To enable automated execution, delete the final line in the script, which waits for user input:
   ```bash
   read -n1 -r -p "Press any key to continue..." key
   ```

These modifications will prepare the script for seamless automated execution.

## Running Your Simulation
Here is the code required to run a DualSPHysics simulation using the Inductiva API:

```python
"""DualSPHysics example."""
import inductiva

# Allocate cloud machine on Google Cloud Platform
cloud_machine = inductiva.resources.MachineGroup(
    provider="GCP",
    machine_type="g2-standard-16",
    spot=True,
    data_disk_gb=200)

# Initialize the Simulator
dualsphysics = inductiva.simulators.DualSPHysics( \
    version="5.4.1")

# Run simulation
task = dualsphysics.run( \
    input_dir=input_dir,
    shell_script="xCaseDambreak3D_FSI_linux64_GPU.sh",
    # Convert VTK files to OBJ files for visualization
    vtk_to_obj=True,
    vtk_to_obj_vtk_dir="CaseDambreak3D_FSI_out/particles/",
    vtk_to_obj_vtk_prefix="PartFluid_",
    vtk_to_obj_particle_radius=0.002,
    vtk_to_obj_smoothing_length=2,
    vtk_to_obj_cube_size=1,
    on=cloud_machine)

# Wait for the simulation to finish and download the results
task.wait()
task.download_outputs()

cloud_machine.terminate()

task.print_summary()
```

This simulation runs in spot mode on a `g2-standard-16` machine, featuring 16 virtual CPUs, 1 NVIDIA L4 GPU, and a 200 GB data disk.

> **Note**: `spot` machines are available at substantial discounts, but your simulation job may be preempted if
> the Cloud provider reclaims the spot machine.

For visualization purposes, this script also optionally converts the simulation’s raw particle data 
(stored as .vtk files) into mesh format (.obj), making it compatible with common visualization tools. 
Check out this [tutorial](https://inductiva.ai/guides/dualsphysics/convert-to-obj) for more details.

When the simulation is complete, we terminate the machine, download the results and print a summary of the simulation as shown below.

```
Task status: Success

Timeline:
	Waiting for Input         at 26/05, 16:53:23      1.581 s
	In Queue                  at 26/05, 16:53:25      67.698 s
	Preparing to Compute      at 26/05, 16:54:32      9.148 s
	In Progress               at 26/05, 16:54:41      6861.774 s
		├> 5862.434 s      bash xCaseDambreak3D_FSI_linux64_GPU.sh
		└> 999.117 s       splashsurf reconstruct CaseDambreak3D_FSI_out/particles//PartFluid_{}.vtk -r=0.002 -l=2 -c=1 -t=0.6 --subdomain-grid=on --mesh-cleanup=on --mesh-smoothing-weights=on --mesh-smoothing-iters=25 --normals=on --normals-smoothing-iters=10 -o CaseDambreak3D_FSI_out/particles//PartFluid__surface{}.obj
	Finalizing                at 26/05, 18:49:03      95.911 s
	Success                   at 26/05, 18:50:39      

Data:
	Size of zipped output:    15.88 GB
	Size of unzipped output:  34.95 GB
	Number of output files:   1230

Estimated computation cost (US$): 0.82 US$
```

As you can see in the "In Progress" line, the part of the timeline that 
represents the actual execution of the simulation, the core computation time 
of this simulation was approximately 1 hour and 37 minutes (5862 seconds).

## Testing Across Different GPUs
As mentioned earlier, the simulation ran on an NVIDIA L4 GPU. To test it on other GPUs, simply change the `machine_type` parameter in the code to `a2-highgpu-1g` or `a3-highgpu-1g`. These machines come equipped with 1 NVIDIA A100 and 1 NVIDIA H100 GPU, respectively.

You also need to specify `zone="europe-west4-b"` when creating the machine, as the A100 and H100 GPUs are not yet available in our default zone.

Here are the results of running the same simulation on these machines:

|  Machine Type  | GPU         |Execution Time          | Estimated Cost |
|:--------------:|:-----------:|:----------------------:|:--------------:|
|  g2-standard-16| Nvidia L4   | 1 hour 36 minutes      | 0.66 US$    |
|  a2-highgpu-1g | Nvidia A100 | 1 hour 2 minutes       | 0.61 US$    |
|  a3-highgpu-1g | Nvidia H100 | 34 minutes 33 seconds  | 1.55 US$    |

> **Note**: The times and costs listed above refer only to running the simulation and do not include the time required to convert VTK files to OBJ format. This conversion step is optional.

Check out our [ParaView](https://inductiva.ai/guides/dualsphysics/paraview-for-visualization) and [Blender](https://inductiva.ai/guides/dualsphysics/blender-for-visualization) tutorials to learn how to visualize your simulation results.
