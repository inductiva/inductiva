# XBeach Visualizations with ParaView
As an example, consider the []() tutorial. -> se já não vier de um tutorial já feito, contextualizar e mostrar como fazer como fizemos aqui:


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

This simulation runs in spot mode on a `g2-standard-16` machine, featuring 16 virtual CPUs, 1 NVIDIA L4 GPU, 
and a 200 GB data disk.

> **Note**: Setting `spot=True` enables the use of spot machines, which are available at substantial discounts. 
> However, your simulation may be interrupted if the cloud provider reclaims the machine.

For visualization purposes, this script also optionally converts the simulation’s raw particle data 
(stored as .vtk files) into mesh format (.obj), making it compatible with common visualization tools. 
Check out this [tutorial](https://inductiva.ai/guides/dualsphysics/convert-to-obj) for more details.

By the end of that tutorial, your simulation will have completed, and the results will be downloaded to the `inductiva_output` folder on your local machine.

This guide will focus on creating a ParaView visualization from those simulation results.

## Visualizing the Results with ParaView
Visualizing your simulation with ParaView is simple and straightforward.

First, open ParaView and go to the menu `File` > `Open...`. Navigate to your
simulation results folder, then to `CaseDambreak3D_FSI_out/particles`, and select the three Groups named `PartFluid_..vtk`, `PartGate_..vtk` and `PartStructure_..vtk`. -> dar update, isto está para o dsph

![File -> Open](./_static/file-open.png)
<p align="center"><em>Figure 1: File -> Open</em></p>

![Selecting the files](./_static/select-files.png) -> dar update, isto está para o dsph
<p align="center"><em>Figure 2: Selecting the files</em></p>

Once all files are loaded, make them visible by clicking the **eye** icon in the **Pipeline Browser** 
on the left side of the screen.

![Make files visible](./_static/eye.png) -> dar update, isto está para o dsph
<p align="center"><em>Figure 3: Make files visible</em></p>

Next, position your camera by clicking the `set view direction +Y` button in the toolbar.

![Move the camera to the correct position](./_static/camera.png)
<p align="center"><em>Figure 4: Move the camera to the correct position</em></p>

Now you can press the **Play** button in the toolbar to watch your simulation run in real time.

![Simulation running](./_static/sim.png) -> dar update, isto está para o dsph
<p align="center"><em>Figure 5: Simulation running</em></p>

## Choosing What Data to Visualize -> dar update, isto está para o dsph
In the previous section, we visualized the particles using ParaView’s default settings. A key part of analyzing your simulation is choosing which data to visualize. For example, DualSPHysics allows you to visualize particle velocity, depending on what data was saved during the simulation.

To do this, select `PartFluid_0000.vtk*` in the **Pipeline Browser** and change the dropdown menu above from `idp` to `Vel`.

![Changing idp to Vel](./_static/pick_vel.png)
<p align="center"><em>Figure 6: Changing idp to Vel</em></p>

![Simulation running with velocity visible](./_static/sim_vel.png)
<p align="center"><em>Figure 7: Simulation running with velocity visible</em></p>

You can save your animation by going to **File > Save Animation...** and choosing your preferred format.