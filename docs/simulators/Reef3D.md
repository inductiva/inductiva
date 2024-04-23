# Reef3D

REEF3D is an open-source hydrodynamics framework with a focus on coastal, marine 
and hydraulic engineering flows. Tailor-made multiphysics solvers are available 
for a range of relevant problems (e.g. sediment transport or floating body dynamics). 
The modular programming approach allows the framework to incorporate a range of 
different flow solvers which together represent all relevant length scales.  
Depending on the wave or flow conditions, the following optimized hydrodynamic 
modules are available:

- **REEF3D::CFD** solves the Navier-Stokes equations in three dimensions. For 
near-field simulations with a complex free surface pattern,  it uses a two-phase 
flow approach with the level set method for interface capturing.
- **REEF3D::FNPF** is a three-dimensional fully nonlinear potential flow solver. 
It is massively parallelized and can be used to create large-scale phase-resolved 
sea states at all water depths.
- **REEF3D::SFLOW** is a depth-averaged model, solving the non-hydrostatic shallow 
water equations ideal for near-shore hydrodynamics and river flow.

## Running a simulation

Reef3D in Inductiva API executes two sequential steps: 
- the meshing with **DiveMESH**;
- the simulation with **Reef3D**. 

Each step is configured with input files, `control.txt` and `ctrl.txt`, respectively. 
Other files may be used to inform the simulator about the grid, geographical data 
or wave information. Reef3D has strict naming policies for each file and we recommend 
users to follow [their guidelines](https://reef3d.wordpress.com/user-guide/). 

To run the simulation users pass a folder with the above files and all others 
required to run the simulation. The folder is then uploaded to Inductiva API and 
the simulation is executed. The output files are then downloaded to the user's 
machine. 

The parallelization of the simulation is handled automatically by Inductiva API, 
based on the number of cores available in the machine.

**General Arguments:**
- `on`: set the machines where the simulations will run. Check [here](../explore_api/computational-infrastructure.md) 
for further detail. If not selected the simulations will be picked-up by a default 
pool shared by everyone.
- `storage_dir`: set the directory where the output files will be stored in the 
cloud. If not selected the output files will be stored in a folder named with the 
task id of the simulation.

For further information on handling the task of the simulation see [here](../explore_api/tasks.md).

## Example

Here, we follow the tutorial with
[regular wave propagation](https://github.com/REEF3D/REEF3D/tree/ed0c8d7a6110892706357f72e0404bd63034efa5/Tutorials/REEF3D_FNPF/9_1%20Regular%20Wave%20Propagation)
from Reef3D repository.

```python
import inductiva

input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "reef3d-input-example.zip", unzip=True
)

reef3d = inductiva.simulators.REEF3D()

task = reef3d.run(input_dir=input_dir)

task.wait()
task.download_outputs()
```
