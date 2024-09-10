# Configuring Simulators

Open-source simulation software lies at the heart of the Indutiva API.
From its inception, the main goal of the API has been to empower scientists
and engineers to take full advantage of state-of-the-art open-source
simulation packages, in a way that minimizes installation issues and leverages
the vast computational power that is currently available in the Cloud. 

At a high level, a significant amount of what we do at Inductiva is to wrap existing 
open-source software packages around a few layers that enable them to execute on
a wide range of virtual machines available on the cloud, and allow simple configuration 
via Python scripting. We wrap such simulation packages in a way that allows us to treat 
them as more abstract computational loads, that have inputs and produce outputs, and 
our job is “merely” just that of passing data around (See [Storage and Data Flow](./data_flow.md)) 
and assigning the simulation tasks to the appropriate computational resource See 
[Resource Allocation Options](./shared_dedicated_resources.md).

Obviously, there is a lot more happening under the hood. For starters, how do we deal 
with the fact that not all simulation software packages work in the same way and, 
therefore, having a fully general formulation for a simulation task is not trivial? 

(the-simple-cases)=
## The simple cases

Some simulation packages offer a single executable that takes as input a single 
file containing the full description of the simulation to be run. This is the simplest 
case because all the API needs to do is send that simulation configuration file to a 
remote VM, and then invoke the simulation executable on that VM, wait for the simulation to 
run, and then forward the results back to the user's remote storage, where they can 
later be retrieved. Of course, installing the simulation package itself may be 
complicated because it may require fulfilling many dependencies, and we deal with all 
of that, such that users don’t have to worry about it. But, from the perspective of the 
orchestration of the simulation, this is quite straightforward.

Our API deals with these cases in a very simple way. In fact, no matter how many
files are required as input for the simulation, we assume that they are all stored
in a folder in the user’s local machine. The entire contents of the folder will be
sent via the API to some remote computational resource. 

Additionally, the executable of the simulator typically expects a “main” simulation configuration
file, from which all other required files (e.g., a 3D model, a bathymetry, etc) are
referred. Below we present the case of two simulators that follow this pattern:

**SWASH:** The main configuration file is a `.sws` file.
```python
# Example of how to run the SWASH simulator
import inductiva

# Set simulation input directory
input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "swash-input-example.zip", unzip=True)

# Instantiate the simulator
swash_simulator = inductiva.simulators.SWASH()

# Input directory contains the .sws config file, a bathymetry file and other files.
task = swash_simulator.run(input_dir=input_dir,
                           sim_config_filename="input.sws")
```

**SPlisHSPlasH:** The main configuration file is a `.json` file.
```python
# Example of how to run the SplishSplash simulator
import inductiva

# Set simulation input directory
input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "splishsplash-input-example.zip", unzip=True)

# Instantiate the simulator
splishsplash_simulator = inductiva.simulators.SplishSplash()

# Input directory contains the .json config file and a .obj file for the domain.
task = splishsplash_simulator.run(input_dir=input_dir,
                                  sim_config_filename="config.json")
```

As you can see, besides the input directory we pass one additional
parameter to the `run()` method: `sim_config_filename`. This refers
to the main configuration file that the simulator executable expects
and for which there is no standard name is expected.

(a-slightly-more-complex-case)=
## A slightly more complex case

Some simulators require running more than one executable to perform a simulation, 
but they are always the same and applied in sequence. In these cases, we 
automatically run both executables together giving the appearance that what is 
being run is a single one. 

An example of this is the REEF3D simulator where the meshing and the simulation
step are run in sequence as a response to a single API command. Below is the
call to the REEF3D simulator that only requires passing the input directory
with the configuration files for both commands.

**REEF3D**

```python
# Example of how to run the REEF3D simulator
import inductiva

# Set simulation input directory
input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "reef3d-input-example.zip", unzip=True
)

# Instantiate the simulator
reef3d_simulator = inductiva.simulators.REEF3D()

# The files for the simulation are in the input directory.
task = reef3d_simulator.run(input_dir=input_dir)
```

In this specific case, REEF3D uses a pre-defined standard for the naming of the
configuration files used by each executable. So, as you can see above, there is
no requirement to pass the `sim_config_filename` parameter. All REEF3D needs is a
pointer to the folder containing all the assets required for the simulation.

(running-long-simulation-pipelines)=
## Running long simulation pipelines

In other simulation packages, a single simulation is more configurable and different
executables can be used to run it. For these cases, users will be able to select  
and configure the executables they wish to run by setting _a priori_ a list of 
commands, where each command describes the executable to be run and the respective
flags to be used. The commands passed are exactly the ones used locally to
run your simulations. 

Further, all the configuration files required by each executable are packed in
a single input directory that is passed to the simulator. The API then uses the 
list of commands and the input directory to run the respective executables in sequence.

Below, an example with OpenFOAM is given.

**OpenFOAM**
```python
# Example of how to run the OpenFOAM simulator
import inductiva

# Set simulation input directory
input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "openfoam-input-example.zip", unzip=True)

# Instantiate the simulator
openfoam_simulator = inductiva.simulators.OpenFOAM()

# Commands to be executed on OpenFOAM
commands = [
    "runApplication surfaceFeatures",
    "runApplication blockMesh",
    "runApplication decomposePar -copyZero",
    "runParallel snappyHexMesh -overwrite",
    "runParallel potentialFoam",
    "runParallel simpleFoam",
    "runApplication reconstructParMesh -constant",
    "runApplication reconstructPar -latestTime"
]

# Run the simulation with the given input directory and commands
task = openfoam_simulator.run(input_dir=input_dir, commands=commands)
```

For this case, the `commands` follow the usual approach used by OpenFOAM
with the `runApplication` and `runParallel` prefix, before stating the
executable, to indicate if the steps are to be run in parallel or not.
All the input files required by each command are set in the input directory.
To run the simulation, all `OpenFOAM` needs is a pointer to the input directory
and the commands.

## What to Read Next

Explore the [available open source simulators](https://docs.inductiva.ai/en/latest/simulators/overview.html)
built into our API.
