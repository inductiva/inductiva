# DualSPHysics simulator

DualSPHysics is a Smoothed-Particle Hydrodynamics (SPH) simulator. The simulator is usually configured by a single file with the extension `.xml`. This file contains all the information about the simulation, including the geometry, the physical properties of the fluids, the boundary conditions, the numerical parameters, and the output files. Sometimes the configuration can also use extra geometry files. 

To run a DualSPHysics simulation you will need to set the commands. In general, there are two main commands that are required to run a simulation: `gencase` and `dualsphysics`. The `gencase` command is used to generate the case files that will be used by the `dualsphysics` command to run the simulation. Thereafter, there are other commands that allow you to post-process the results. Examples of this are:
- `partvtk`: to generate VTK files with the particle trajectories;
- `isosurface`: to generate VTK files with the isosurfaces of the fluid;
- `measuretool`: to generate CSV files with measurements of the fluid properties.

For an extensive list of commands please check the DualSPHysics [documentation](https://dual.sphysics.org/). For the API, you can passed them in lower case and we will handle the rest for you!

## Example

In this example, we run a classical CFD case of a flow over a cylinder. 

```python
import inductiva

# Download the configuration files into a folder
input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "dualsphysics-input-example.zip", unzip=True
)

commands = [{"cmd": "gencase config flow_cylinder -save:all", "prompts": []},
{"cmd": "dualsphysics flow_cylinder flow_cylinder -dirdataout data -svres", "prompts": []},
{"cmd": "partvtk -dirin flow_cylinder/data -savevtk flow_cylinder/PartFluid -onlytype:-all,+fluid", "prompts": []}]

# Initialize the Simulator
dualsphysics = inductiva.simulators.DualSPHysics()

# Run simulation with config files in the input directory
task = dualsphysics.run(input_dir=input_dir,
                        commands=commands)

task.wait()
task.download_outputs()

```