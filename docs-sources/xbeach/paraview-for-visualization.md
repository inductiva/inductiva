# XBeach Visualizations with ParaView
This tutorial guides you through creating a visualization in ParaView using simulation results from XBeach.

As an example, we will use the `DELILAH` field experiment, which is also featured in the [official XBeach documentation](https://xbeach.readthedocs.io/en/stable/examples.html#field-experiment-delilah).

## Prerequisites
Download the necessary input files from this [link](https://svn.oss.deltares.nl/repos/xbeach/skillbed/input/Delilah_199010131000/). 

### Visualization requirements — DELILAH case

In order to automatically generate the XBeach visualization at the end of your DELILAH simulation ensure that the following variables exists in the `params_original.txt` input file, or modify them accordingly:

- `outputformat = netcdf`: Ensure that the xbeach is exporting the results in the `netcdf` format
- `nglobalvar = 3` followed by `H`, `Zs`, and `Zb`. These three global variables gives the visualization script everything it needs to rebuild the water surface (H, Zs) and the seabed (Zb) in 3D.
- `single_dir = 0`: – required for this particular example.
- `tintg = 5`: tells XBeach to write those global variables every few time steps; the low value generates enough frames for a smooth animation without creating huge files.



## Running the DELILAH case
Below is a script to run the DELILAH case using the Inductiva API. It is configured to export simulation results in a .vtk format, via the `export_vtk` flag, which is compatible with ParaView for visualization.

```python
import inductiva

# Allocate cloud machine on Google Cloud Platform
cloud_machine = inductiva.resources.MachineGroup(
    provider="GCP",
    machine_type="c2d-highcpu-16",
    spot=True)

# Set simulation input directory
input_dir = "Path/to/Delilah_199010131000"

# Initialize the Simulator
xbeach = inductiva.simulators.XBeach()

# Run simulation
task = xbeach.run(
    input_dir=input_dir,
    sim_config_filename="params_original.txt",
    on=cloud_machine,
    project="xbeach",
    export_vtk=True, # Flag to control whether to generate the visualization 
)

# Wait for the simulation to finish and download the results
task.wait()
task.download_outputs()

cloud_machine.terminate()

task.print_summary()
```

At the end of the run, the simulation will be complete, and the results will be downloaded to the `inductiva_output` folder on your local machine.


## Visualizing the Results with ParaView

Set the `export_vtk` flag to `True` then a `VTK` folder will be created in the simulation output directory. Inside of this folder there should be two vtk groups: 

1. seabed.vtk -> single, static mesh of the seabed.
2. wave_..vtk -> a group of .vtk files which contains water-surface data at every `tintg` interval.

Visualizing your simulation with ParaView is simple and straightforward.

First, open ParaView and go to the menu `File` > `Open...`. Navigate to your
simulation results folder, then to `VTK/`, and select both Groups named `seabed.vtk`, and `wave_..vtk`.

![File -> Open](./_static/file-open.png)
<p align="center"><em>Figure 1: File -> Open</em></p>

![Selecting the files](./_static/select-files.png)
<p align="center"><em>Figure 2: Selecting the files</em></p>

Once all files are loaded, make them visible by clicking the **eye** icon in the **Pipeline Browser** 
on the left side of the screen.

![Make files visible](./_static/eye.png)
<p align="center"><em>Figure 3: Make files visible</em></p>

Next, position your camera by clicking the `set view direction +X` button in the toolbar.

![Move the camera to the correct position](./_static/camera.png)
<p align="center"><em>Figure 4: Move the camera to the correct position</em></p>

Now you can press the **Play** button in the toolbar to watch your simulation run in real time.

![Simulation running](./_static/sim.gif)
<p align="center"><em>Figure 5: Simulation running</em></p>

You can save your animation by going to **File > Save Animation...** and choosing your preferred format.

That's it! You’ve completed the full workflow: from running an XBeach simulation using the DELILAH case to visualizing the results 
in ParaView.