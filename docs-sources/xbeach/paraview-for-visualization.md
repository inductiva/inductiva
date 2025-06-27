# XBeach Visualizations with ParaView
This tutorial guides you through creating a visualization in ParaView using simulation results from XBeach.

As an example, we will use the `DELILAH` field experiment, which is also featured in the [official XBeach documentation](https://xbeach.readthedocs.io/en/stable/examples.html#field-experiment-delilah).

## Prerequisites
Download the necessary input files from this [link](https://svn.oss.deltares.nl/repos/xbeach/skillbed/input/Delilah_199010131000/). 

Alguma versão do Paraview que devamos recomendar?

## Running the DELILAH case
Below is a script to run the DELILAH case using the Inductiva API. It is configured to export simulation results in a .vtk format, which is compatible with ParaView for visualization.

```python
import inductiva

# Allocate cloud machine on Google Cloud Platform
cloud_machine = inductiva.resources.MachineGroup(
    provider="GCP",
    machine_type="c2d-highcpu-16",
    spot=True)

# Initialize the Simulator
dualsphysics = inductiva.simulators.DualSPHysics( \
    version="1.24")


# Run simulation
COLOCAR TASK + FLAG export_vtk = True

# Wait for the simulation to finish and download the results
task.wait()
task.download_outputs()

cloud_machine.terminate()

task.print_summary()
```

At the end of the run, the simulation will be complete, and the results will be downloaded to the `inductiva_output` folder on your local machine.

You’re now ready to start creating your visualization!

## Visualizing the Results with ParaView
Visualizing your simulation with ParaView is simple and straightforward.

First, open ParaView and go to the menu `File` > `Open...`. Navigate to your
simulation results folder, then to `CaseDambreak3D_FSI_out/particles`, and select the three Groups named `PartFluid_..vtk`, `PartGate_..vtk` and `PartStructure_..vtk`. -> DAR UPDATE

![File -> Open](./_static/file-open.png)
<p align="center"><em>Figure 1: File -> Open</em></p>

![Selecting the files](./_static/select-files.png) -> DAR UPDATE, isto está para o DSPH
<p align="center"><em>Figure 2: Selecting the files</em></p>

Once all files are loaded, make them visible by clicking the **eye** icon in the **Pipeline Browser** 
on the left side of the screen.

![Make files visible](./_static/eye.png) -> DAR UPDATE
<p align="center"><em>Figure 3: Make files visible</em></p>

Next, position your camera by clicking the `set view direction +Y` button in the toolbar.

![Move the camera to the correct position](./_static/camera.png)
<p align="center"><em>Figure 4: Move the camera to the correct position</em></p>

Now you can press the **Play** button in the toolbar to watch your simulation run in real time.

![Simulation running](./_static/sim.png) -> DAR UPDATE
<p align="center"><em>Figure 5: Simulation running</em></p>

## Choosing What Data to Visualize -> DAR UPDATE 
In the previous section, we visualized the particles using ParaView’s default settings. A key part of analyzing your simulation is choosing which data to visualize. For example, DualSPHysics allows you to visualize particle velocity, depending on what data was saved during the simulation.

To do this, select `PartFluid_0000.vtk*` in the **Pipeline Browser** and change the dropdown menu above from `idp` to `Vel`.

![Changing idp to Vel](./_static/pick_vel.png)
<p align="center"><em>Figure 6: Changing idp to Vel</em></p>

![Simulation running with velocity visible](./_static/sim_vel.png)
<p align="center"><em>Figure 7: Simulation running with velocity visible</em></p>

You can save your animation by going to **File > Save Animation...** and choosing your preferred format.

That's it! You’ve completed the full workflow: from running an XBeach simulation using the DELILAH case to visualizing the results 
in ParaView.