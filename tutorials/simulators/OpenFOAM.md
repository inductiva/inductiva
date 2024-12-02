In this guide, we will walk you through setting up and running OpenFOAM 
simulations using the Inductiva API.

We will cover:

- Configuring OpenFOAM simulations with the appropriate input directories.
- Example codes to help you get started with simulations.
- An advanced example for running MB9 Micro-benchmark by ExaFOAM.

# OpenFOAM

OpenFOAM is a Finite Volume method for CFD simulations with a wide range of 
applications across several areas of engineering and science. It offers 
a broad set of features for everything from **complex fluid flows** (including 
chemical reactions, turbulence, and heat transfer) to **solid dynamics** and 
**electromagnetics**.

There are two main open-source distributions of OpenFOAM: one developed by the
[OpenFOAM foundation](https://openfoam.org/) and another by the
[ESI Group](https://www.openfoam.com/). The Inductiva API supports both,
and you can select your preferred distribution by setting the `distribution` parameter
when initializing the simulator. *By default, it uses the OpenFOAM Foundation version.*

A single simulation via the Inductiva API follows several steps in OpenFOAM, 
including domain partitioning, meshing, solving, and post-processing. To 
configure a simulation for OpenFOAM, you’ll need a set of configuration 
files organized into three folders:
- `time`: Contains files for particular fields, like
initial values and boundary conditions that you must specify. For example, 
initial conditions at  t=0  are stored in the `0` directory.
- `constant`: Contains files that describe the objects in the simulation and the 
physical properties being modeled.
- `system`: Contains files that describe the simulation, including solvers, 
numerical parameters, and output files. 

It must contain at least 3 files: 
`controlDict`: Run control parameters (start/end time, time step, and parameters 
for data output)
`fvSchemes`: Selection of discretization schemes used during the solution.
`fvSolution`: Equation solvers, tolerances and other algorithm controls.

All of these folders should be placed inside an **input directory**. To run a 
simulation, you need to configure a list of dictionaries that specify the 
commands to be executed on the backend. 

Below is an example of how to run the [motorbike tutorial](https://github.com/OpenFOAM/OpenFOAM-8/tree/master/tutorials/incompressible/simpleFoam/motorBike) 
from OpenFOAM, demonstrating how this process works in practice.

## Example Code - OpenFOAM Foundation Distribution

In this example, we demonstrate how to run the [motorbike tutorial](https://github.com/OpenFOAM/OpenFOAM-8/tree/master/tutorials/incompressible/simpleFoam/motorBike) 
tutorial using the OpenFOAM Foundation distribution.

```{literalinclude} ../../examples/openfoam_foundation/openfoam_foundation.py
:language: python
```

## Example Code - ESI Distribution

To run the sample simulation above, simply download the
`openfoam-esi-input-example.zip` file, select the correct distribution by
using `inductiva.simulators.OpenFOAM(distribution="esi")`, and run the respetive
`Allrun` file.

## Set the commands manually

If you decide to set your commands manually, you can and here are a two examples:

```python
# Commands to run using single machine
commands_single_machine = [
    "runApplication surfaceFeatures",
    "runApplication blockMesh",
    "runApplication decomposePar -copyZero",
    "runParallel snappyHexMesh -overwrite",
    "runParallel potentialFoam",
    "runParallel simpleFoam",
    "runApplication reconstructParMesh -constant",
    "runApplication reconstructPar -latestTime"
]
```

```python
# Commands to run using cluster
config = inductiva.commands.MPIConfig("4.1.6",np=4,use_hwthread_cpus=False)
commands_cluster = [
    "runApplication surfaceFeatures",
    "runApplication blockMesh",
    "runApplication decomposePar -copyZero",
    inductiva.commands.Command("snappyHexMesh -overwrite -parallel", mpi_config=config),
    inductiva.commands.Command("potentialFoam -parallel", mpi_config=config),
    inductiva.commands.Command("simpleFoam -parallel", mpi_config=config),
    "runApplication reconstructParMesh -constant",
    "runApplication reconstructPar -latestTime"
]
```

In the first example, the MPI configuration is managed automatically by the
OpenFOAM simulator.

In the second example, you configure the MPI settings manually. This approach
allows you to select the MPI version, specify the number of processes, and
enable or disable hyperthreading based on your requirements.  

One advantage of the second approach is that by explicitly defining the MPI
configuration through our API, you can set up MPI clusters and leverage multiple
machines for running your simulations.  

For more details on commands and MPI configuration, refer to the
[Custom Docker Images](https://tutorials.inductiva.ai/simulators/CustomImage.html#command-and-mpiconfig)
documentation.

## Advanced Tutorial: Running the MB9 Micro-benchmark from ExaFOAM

This guide walks you through running a complex OpenFOAM simulation using the
**MB9 micro-benchmark** from [ExaFOAM](https://exafoam.eu/benchmarks/). This
benchmark simulates a high-lift aircraft configuration, ideal for studying
near-wall turbulence using **wall-modeled Large Eddy Simulation (WMLES)**.

### Objective

We'll run this simulation on:
1. **Single 360 vCPU Machine**.
2. **MPI Cluster** using two 360 vCPU machines to improve performance.

### Prerequisites

1. **Download Input Files**: Get the input files from the
[repository](https://develop.openfoam.com/committees/hpc/-/tree/develop/compressible/rhoPimpleFoam/LES/highLiftConfiguration)
and place them in a folder named `highLiftConfiguration`.

   **Directory Structure**:
   ```bash
   ls -lasgo highLiftConfiguration
   total 104
   0 drwxrwxr-x@ 14     448 Sep 23 11:44 .
   0 drwx------@ 19     608 Sep 23 11:49 ..
   0 drwxrwxr-x@ 12     384 Sep 20 09:43 0.orig
   8 -rwxr-xr-x@  1     626 Jun 13 10:09 Allclean
   16 -rwxr-xr-x@  1    6998 Jun 13 10:09 Allrun
   8 -rw-rw-r--@  1     991 Jun 13 10:09 COPYING
   48 -rw-rw-r--@  1   21547 Jun 13 10:09 README.md
   0 -rw-rw-r--@  1       0 Jun 13 10:09 case.foam
   0 drwxrwxr-x@  6     192 Sep 20 09:43 constant
   0 drwxrwxr-x@ 13     416 Jun 13 10:09 figures
   0 drwxrwxr-x@ 28     896 Sep 20 09:43 system
   24 -rw-rw-r--@  1   11399 Jun 13 10:09 thumbnail.png
   ```

### Overview

Here’s the code you'll be working on as we progress through the tutorial. Don’t
worry if it doesn’t all make sense right now; everything will become clearer
in the upcoming steps.

```python
import inductiva

machine_group = inductiva.resources.MachineGroup(
            machine_type="c3d-highcpu-360",
            spot=True)
machine_group.start()

input_dir = "/path/to/highLiftConfiguration"

#Choose your simulator
openfoam = inductiva.simulators.OpenFOAM(distribution="esi")

task = openfoam.run(
            input_dir=input_dir,
            commands=["bash ./Allrun"],
            n_vcpus=180,
            use_hwthread=True,
            on=machine_group)

task.wait()
task.download_outputs()
machine_group.terminate()

task.print_summary()

```

### Step 1: Adjust Simulation Parameters

For a faster simulation, modify the following parameters in the case definition
file (`system/include/caseDefinition`):

- **Time Step (`dt`)**: Set to 0.00002.
- **Start Time (`initTime`)**: 0.10
- **End Time (`finalTime`)**: 0.30

### Step 2: Running the Simulation

#### a. Configure and Start Machine

1. **Pick your machine**:
    ```python
    import inductiva
    machine_group = inductiva.resources.MachineGroup(machine_type="c3d-highcpu-360", spot=True)
    ```
    **Note**: `spot` machines are a lot cheaper but can be terminated by the
    provider if needed.

2. **Start your machine**
    ```python
    machine_group.start()
    ```

#### b. Simulation inputs
1. **Specify Simulation Directory**:
Let's start by defining a variable that points to the `highLiftConfiguration`
folder where all your simulation files are located.

   ```python
   input_dir = "/path/to/highLiftConfiguration"
   ```
2. **Read Commands**:
Now, to run this simulation you need only to simply run the `Allrun` file.
   ```python
   commands = ["bash ./Allrun"]
   ```

#### c. Run your simulation

1. **Run the simulation**:
We now have all we need to run our simulation.
   ```python
   #Choose your simulator
   openfoam = inductiva.simulators.OpenFOAM(distribution="esi")

   task = openfoam.run(
               input_dir=input_dir,
               commands=["bash ./Allrun"],
               on=machine_group)
   ```

2. **Wait and Download Outputs**:
That is it. Our simulation is now running on the cloud. We can `wait` for the
simulation to be over, or we can turn our computer off go for a coffe (☕️).
   ```python
   task.wait()
   task.download_outputs()
   ```

3. **Terminate Machine**:
Once our simulation is over we can/should terminate our machine to save on costs.
If you forget, dont worry we got your back. By default, a machine will be
automaticly terminated if no simulation runs on it for 30 minutes.

   ```python
   machine_group.terminate()
   ```

4. **Check your simulation summary**:
Now that our simulation has finished we can print a summary of said simulation.
This includes information about the execution times, outputs generated and
much more.
    ```python
    task.print_summary()
    ```

### Step 4: Enhancing Performance with MPI Cluster

Since we now want to scale our simulation into multiple machines we can no longer
let Openfoam take care of our MPI configuration. We need to do it manually to
let our API know how to handle MPI.

The first thing we need to do is move from running `Allrun` into setting the
commands manually with the correct MPI configuration.

1. **MPI Configuration**:

    The first thing to do is defining the MPI configuration. This is done by:

    ```python
    config = inductiva.commands.MPIConfig("4.1.6",np=180,use_hwthread_cpus=False)
    ```

    Apart from that, we need to edit the `nCores` to 180 in the `system/include/caseDefinition`
    file.

2. **Commands**:

    Now we need to define the commands to run on the cluster. This is done by
    "translating" the `Allrun` file into a list of commands. Like so:

    ```python
    commands = [
        'cp system/controlDict.SHM system/controlDict',
        'cp system/fvSchemes.SHM system/fvSchemes',
        'cp system/fvSolution.SHM system/fvSolution',
        'runApplication blockMesh',
        'runApplication blockMesh',
        'runApplication snappyHexMesh -dict system/snappyHexMeshDict.refineblockMesh -overwrite',
        'mv ./0/cellLevel ./constant/polyMesh/',
        'mv ./0/pointLevel ./constant/polyMesh/',
        'rm -r 0',
        'runApplication checkMesh',
        'runApplication extrudeMesh -dict system/extrudeMeshDict.refineblockMesh',
        'rm constant/polyMesh/cellLevel',
        'rm constant/polyMesh/cellZones',
        'rm constant/polyMesh/pointLevel',
        'rm constant/polyMesh/pointZones',
        'rm constant/polyMesh/faceZones',
        'rm constant/polyMesh/level0Edge',
        'rm constant/polyMesh/surfaceIndex',
        'sed -i s/"ff_zMin"/"ff_SAVE"/g constant/polyMesh/boundary',
        'sed -i s/"ff_zMax"/"ff_zMin"/g constant/polyMesh/boundary',
        'sed -i s/"ff_SAVE"/"ff_zMax"/g constant/polyMesh/boundary',
        'runApplication checkMesh',
        'runApplication surfaceFeatureExtract',
        'runApplication decomposePar -decomposeParDict system/decomposeParDict.SHM',
        inductiva.commands.Command('snappyHexMesh -decomposeParDict system/decomposeParDict.SHM -dict system/snappyHexMeshDict -overwrite',mpi_config=config),
        inductiva.commands.Command('checkMesh -decomposeParDict system/decomposeParDict.SHM -latestTime -meshQuality',mpi_config=config),
        'runApplication reconstructParMesh -mergeTol 1e-08 -constant -latestTime',
        'cp -r constant/polyMesh constant/polyMesh.origHalf',
        'runApplication mirrorMesh -overwrite',
        'runApplication topoSet -dict system/topoSetDict.faces.ff_zMin',
        'runApplication createPatch -overwrite -dict system/createPatchDict.mirrorMesh',
        'runApplication changeDictionary -constant -dict system/changeDictionaryDict.cyclicPatches -enableFunctionEntries',
        'runApplication topoSet -dict system/topoSetDict.faces.cyclic',
        'runApplication checkMesh -constant',
        'rm -r processor*',
        'rm -r constant/polyMesh.origHalf',
        'cp -r 0.orig 0',
        'cp system/controlDict.SHM system/controlDict',
        'cp system/fvSchemes.SHM system/fvSchemes',
        'cp system/fvSolution.SHM system/fvSolution',
        'runApplication decomposePar',
        inductiva.commands.Command('renumberMesh -overwrite',mpi_config=config),
        'cp system/controlDict.SRS.init system/controlDict',
        'cp system/fvSolution.SRS system/fvSolution',
        'cp system/fvSchemes.SRS system/fvSchemes',
        inductiva.commands.Command('applyBoundaryLayer -ybl 0.1',mpi_config=config),
        'cp system/controlDict.SRS.init system/controlDict',
        'cp system/fvSolution.SRS system/fvSolution',
        'cp system/fvSchemes.SRS system/fvSchemes',
        inductiva.commands.Command('rhoPimpleFoam',mpi_config=config),
        'cp system/controlDict.SRS.avg system/controlDict',
        'cp system/fvSolution.SRS system/fvSolution',
        'cp system/fvSchemes.SRS system/fvSchemes',
        inductiva.commands.Command('rhoPimpleFoam',mpi_config=config),
        inductiva.commands.Command('postProcess -func sampleDict.surface.SRS -latestTime',mpi_config=config),
    ]
    ```
3. **Start your Cluster**:

    Now we need to start our cluster. This is done by:

    ```python
    mpi_cluster = inductiva.resources.MPICluster(
                    machine_type="c3d-highcpu-180",
                    data_disk_gb=300,
                    num_machines=2)
    mpi_cluster.start()
    ```
    
4. **Run the simulation**:
    All that is left is to run the simulation on the cluster.

    ```python
    task = openfoam.run(
                    input_dir=input_dir,
                    commands=commands,
                    on=mpi_cluster)
    ```

As you can see the process of scalling up (or down) can be done easly by just
picking a new resource. We encorage you to try other machines/configurations.

**Note**: `spot` machines are not supported for MPI clusters.

### Conclusion

Running the simulation on a high-performance machine and scaling it on an MPI
cluster can significantly reduce computation time. Download and analyze your
results locally once complete.

Happy simulations!

## What to read next

Try our [Inductiva API](https://console.inductiva.ai/) 
to streamline your workflows and make the most of cloud resources for 
large-scale simulations.

You may also be interested in reading our blog post,
[The 3D Mesh Resolution Threshold - 5k Points is All You Need!](https://inductiva.ai/blog/article/5k-points-is-all-you-need), 
where we explore just how much you can reduce the level of detail in a 
3D object while still maintaining accurate aerodynamic results in a virtual 
wind tunnel built with OpenFOAM.
