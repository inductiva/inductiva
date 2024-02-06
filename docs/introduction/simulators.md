# Simulators

Open-source simulation software lies at the heart of the Indutiva API. From its
inception, the main goal of the API has been to empower scientists and engineers
to take full advantage of state-of-the-art open-source simulation packages, in a
way that minimizes installation issues and leverages the vast computational power
that is currently available in the Cloud. 

At a high level, a significant amount of what we do at Inductiva is to wrap existing 
open-source software packages around a few layers that enable them to execute on
a wide range of virtual machines available on the cloud, and allow simple configuration 
via Python scripting. We wrap such simulation packages in a way that allows us to treat 
them as more abstract computational loads, that have inputs and produce outputs, and 
our job is “merely” just that of passing data around (See [Storage and Data Flow]()) 
and assigning the simulation tasks to the appropriate computational resource (See 
[Shared and Dedicated Resources]()).

Obviously, there is a lot more happening under the hood. For starters, how do we deal 
with the fact that not all simulation software packages work in the same way and, 
therefore, having a fully general formulation for a simulation task is not trivial? 

## The simple cases

Some simulation packages offer a single binary that takes as input a single 
file containing the full description of the simulation to be run. This is the simplest 
case because all the API needs to do is send that simulation configuration file to a 
remote VM, and then invoke the simulation binary on that VM, wait for the simulation to 
run, and then forward the results back to the user's remote storage, where they can 
later be retrieved. Of course, installing the simulation package itself may be 
complicated because it may require fulfilling many dependencies, and we deal with all 
of that, such that users don’t have to worry about it. But, from the perspective of the 
orchestration of the simulation, this is quite straightforward.

Our API deals with these cases in a very simple way. In fact, no matter how many
files are required as input for the simulation, we assume that they are all stored
in a folder in the user’s local machine. The entire contents of the folder will be
sent via the API to some remote computational resource. For some simulators that is
all that is needed, e.g., REEF3D. 


For others, a specific binary of the simulator expects a “main” simulation configuration
file, from which all other required files (e.g., a 3D model, a bathymetry, etc) are
referred. Below we present two simulators that follow this pattern:

```python
# Example of how to run the SWASH simulator
import inductiva

# Instantiate the simulator
swash_simulator = inductiva.simulators.SWASH()

# SWASH requires a .sws file as the main configuration file.
# Input directory contains the main config file, the bathymetry file and other
# files configuring the domain
task = swash_simulator.run(input_dir="swash-example",
                           sim_config_filename="input.sws")

```

```python
# Example of how to run the SplishSplash simulator
import inductiva

# Instantiate the simulator
splishsplash_simulator = inductiva.simulators.SplishSplash()

# SPlisHSPlasH requires a .json file as the main config file.
task = splishsplash_simulator.run(input_dir="splishsplash-example",
                                  sim_config_filename="config.json")
```

## A slightly more complex case

Some simulators require running more than one command, but they are
always the same set of commands in sequence. In these cases, we automatically run both 
commands giving the appearance that what is being run is a single command. An example of this is the REEF3D simulator where the meshing and the simulation step
are run in sequence as a response to a single API command. Below is the
call to the REEF3D simulator that only requires passing the input directory
with the configuration files for both commands.

```python
# Example of how to run the REEF3D simulator
import inductiva

# Instantiate the simulator
reef3d_simulator = inductiva.simulators.REEF3D()

# REEF3D only needs the input directory to be passed to run. The files for the
# simulation run are all contained within the input directory.
task = reef3d_simulator.run(input_dir="reef3d-example")
```

## Running long simulation pipelines

In other simulation packages, a single simulation can be composed of running
several binaries or commands in sequence. In these cases, to run the simulation
the configuration files used by each command are contained in a single input directory.
This one is passed to the simulator as before, with the addition of a list of
commands to be run. The commands are treated as if the simulation running
locally. An example is given in the next code snippet with OpenFOAM.

```python
# Example of how to run the OpenFOAM simulator
import inductiva

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

# OpenFOAM requires the commands to be executed and the input directory with all 
# configuration files.
task = openfoam.run(input_dir=input_dir, commands=commands)
```

## Available Simulator Packages

Inductiva API has available several open-source simulators ready to use. Users 
who are familiar with the simulators can easily start running simulations with 
their previously prepared simulation configuration files. 

The simulators currently available are:
- [SPlisHSPlasH](./simulators/SPlisHSPlasH.md)
- [DualSPHysics](./simulators/DualSPHysics.md)
- [OpenFOAM](./simulators/OpenFOAM.md)
- [SWASH](./simulators/SWASH.md)
- [XBeach](./simulators/XBeach.md)
- [Reef3D](./simulators/Reef3D.md)
- [GROMACS](./simulators/GROMACS.md)
- [FDS](./simulators/FDS.md)

Check the documentation of each simulator to learn more on how to configure them. 
