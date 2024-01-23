## OpenFOAM simulator

OpenFOAM is a Finite Volume method for CFD simulations with a wide range of 
applications across several areas of engineering and science. OpenFOAM has an 
extensive
range of features to solve anything from complex fluid flows involving chemical 
reactions, turbulence and heat transfer, to solid dynamics and electromagnetics.

There are two main open-source versions of OpenFOAM, one developed by [OpenFOAM foundation](https://openfoam.org/) and another by the [ESI Group](https://www.openfoam.com/). Inductiva API supports 
both versions, and you can choose which one you want to use by setting the `version` 
parameter when initializing the simulator. The default version is the one developed 
by the OpenFOAM foundation.

A single simulation via Inductiva API comprises several steps done via OpenFOAM 
- e.g., partitioning the domain, meshing, solvers and post-processing. Hence, to 
configure a simulation for OpenFOAM the user will need a set of configuration files 
that are organized in three folders:
- `time`: containing individual files of data for particular fields, like initial 
values and boundary conditions that the user must specify. For example, for an 
initial condition at $, the initial conditions will be stored in the directory `0`.
- `constant`: contains files that describe the objects in the simulation and the 
physical properties of the application we are concerned.
- `system`: contains all of the files that describe the simulation, including the 
solvers, the numerical parameters, and the output files. It must contain at least 
3 files: `controlDict` where run control parameters are set including start/end time, 
time step and parameters for data output; `fvSchemes` where discretization schemes 
used in the solution may be selected at run-time; and `fvSolution` where the equation 
solvers, tolerances and other algorithm controls are set for the run.

All of these folders should be inside an input directory. Finally, to run a 
simulation the user needs to configure a list of dictionaries specifying the commands 
they want to execute on the backend. Below, we run the [motorbike tutorial](https://github.com/OpenFOAM/OpenFOAM-8/tree/master/tutorials/incompressible/simpleFoam/motorBike) from OpenFOAM and show 
how this is done in practice.

The commands passed to the simulator follow the structure of OpenFOAM, that is, 
using the prefix `runApplication` the command will execute sequentially and with
 `runParallel` the command will use the maximum number of cores the machine has 
 available. Hence, you don't need to set the specific number of processes, the 
 simulator will do that for you automatically. In particular, the decomposeParDict 
 will be configured automatically and, at the moment, only the scotch decomposition 
 method is available.

### Example 

````python
import inductiva

# Set simulation input directory
input_dir = inductiva.utils.files.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "openfoam-input-example.zip", unzip=True)

# Set the simulation commands
commands = parallel_commands = [
    {"cmd": "runApplication surfaceFeatures", "prompts": []},
    {"cmd": "runApplication blockMesh", "prompts":[]},
    {"cmd": "runApplication decomposePar -copyZero", "prompts":[]},
    {"cmd": "runParallel snappyHexMesh -overwrite", "prompts":[]},
    {"cmd": "runParallel potentialFoam", "prompts":[]},
    {"cmd": "runParallel simpleFoam", "prompts":[]},
    {"cmd": "runApplication reconstructParMesh -constant", "prompts":[]},
    {"cmd": "runApplication reconstructPar -latestTime", "prompts": []}
]

# Initialize the Simulator
openfoam = inductiva.simulators.OpenFOAM(version="foundation")

# Run simulation with config files in the input directory
task = openfoam.run(input_dir=input_dir, commands=commands)

task.wait()
task.download_outputs()
````