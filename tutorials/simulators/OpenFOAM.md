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

The commands passed to the simulator follow OpenFOAM’s structure. Using 
the `runApplication` prefix will execute commands sequentially, while
`runParallel` will use all available CPU cores automatically—no need to 
manually set the number of processes. The **decomposeParDict** is configured 
automatically and currently, only the **scotch decomposition method** is supported.

## Example Code - OpenFOAM Foundation Distribution

In this example, we demonstrate how to run the [motorbike tutorial](https://github.com/OpenFOAM/OpenFOAM-8/tree/master/tutorials/incompressible/simpleFoam/motorBike) 
tutorial using the OpenFOAM Foundation distribution.

```{literalinclude} ../../examples/openfoam_foundation/openfoam_foundation.py
:language: python
```

## Example Code - ESI Distribution

To run the sample simulation above, simply download the
`openfoam-esi-input-example.zip` file, select the correct distribution by
using `inductiva.simulators.OpenFOAM(distribution="esi")`, and replace  
`runApplication surfaceFeatures` with `runApplication surfaceFeatureExtract`.

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

### Step 1: Adjust Simulation Parameters

For a faster simulation, modify the following parameters in the case definition
file (`system/include/caseDefinition`):

- **Time Step (`dt`)**: Set to 0.00002.
- **Start Time (`initTime`)**: 0.10
- **End Time (`finalTime`)**: 0.30

### Step 2: Create Commands File

For now, we only allow the user to send a **list of specific commands** to the
simulator. This commands are `runApplication` and `runParallel`.

Due to this limitation,we need to "translate" the `Allrun` file into a list of
commands that can be sent to the simulator.

Another important point is that the `Allrun` file contains a lot of `cd`
commands, which are not supported by the simulator. As a workaround, we can
replace the `cd` commands with `runApplication` followed by `cd` (the same
applies to `mv` and other commands).

As an example, the `Allrun` file contains the following commands:

```bash
snappyHexMesh -dict system/snappyHexMeshDict.refineblockMesh -overwrite > ${LOGDIR}/log.M02.snappyHexMesh.refineblockMesh 2>&1 || exit 1
mv ./0/cellLevel ./constant/polyMesh/
mv ./0/pointLevel ./constant/polyMesh/

```
Those commands can be translated into the following:

```bash
runApplication snappyHexMesh -dict system/snappyHexMeshDict.refineblockMesh -overwrite
runApplication mv ./0/cellLevel ./constant/polyMesh/
runApplication mv ./0/pointLevel ./constant/polyMesh/
```

Since we only allow a list of commands all `if` statements in the `Allrun` file
should be removed. And a clear line of execution should be defined.

```bash
if condition ; then
    runApplication command1
    runApplication command2
else
    runApplication command3
    runApplication command4
fi
```

Should be translated to:

- If we want to simulate the `if` condition

```bash
runApplication command1
runApplication command2
```

- If we want to simulate the `else` condition

```bash
runApplication command3
runApplication command4
```

In the end, the `Allrun` file should be translated into a list of commands that
should be placed into a file named `commands.txt`.

Said file should contain the following commands:

```bash
runApplication cp system/controlDict.SHM system/controlDict
runApplication cp system/fvSchemes.SHM system/fvSchemes
runApplication cp system/fvSolution.SHM system/fvSolution
runApplication blockMesh
runApplication blockMesh
runApplication snappyHexMesh -dict system/snappyHexMeshDict.refineblockMesh -overwrite
runApplication mv ./0/cellLevel ./constant/polyMesh/
runApplication mv ./0/pointLevel ./constant/polyMesh/
runApplication rm -r 0
runApplication checkMesh
runApplication extrudeMesh -dict system/extrudeMeshDict.refineblockMesh
runApplication rm constant/polyMesh/cellLevel
runApplication rm constant/polyMesh/cellZones
runApplication rm constant/polyMesh/pointLevel
runApplication rm constant/polyMesh/pointZones
runApplication rm constant/polyMesh/faceZones
runApplication rm constant/polyMesh/level0Edge
runApplication rm constant/polyMesh/surfaceIndex
runApplication sed -i s/"ff_zMin"/"ff_SAVE"/g constant/polyMesh/boundary
runApplication sed -i s/"ff_zMax"/"ff_zMin"/g constant/polyMesh/boundary
runApplication sed -i s/"ff_SAVE"/"ff_zMax"/g constant/polyMesh/boundary
runApplication checkMesh
runApplication surfaceFeatureExtract
runApplication decomposePar -decomposeParDict system/decomposeParDict.SHM
runParallel snappyHexMesh -decomposeParDict system/decomposeParDict.SHM -dict system/snappyHexMeshDict -overwrite
runParallel checkMesh -decomposeParDict system/decomposeParDict.SHM -latestTime -meshQuality
runApplication reconstructParMesh -mergeTol 1e-08 -constant -latestTime
runApplication cp -r constant/polyMesh constant/polyMesh.origHalf
runApplication  mirrorMesh -overwrite
runApplication topoSet -dict system/topoSetDict.faces.ff_zMin
runApplication createPatch -overwrite -dict system/createPatchDict.mirrorMesh
runApplication changeDictionary -constant -dict system/changeDictionaryDict.cyclicPatches -enableFunctionEntries
runApplication topoSet -dict system/topoSetDict.faces.cyclic
runApplication checkMesh -constant
runApplication rm -r processor*
runApplication rm -r constant/polyMesh.origHalf
runApplication cp -r 0.orig 0
runApplication cd system
runApplication cp system/controlDict.SHM system/controlDict
runApplication cp system/fvSchemes.SHM system/fvSchemes
runApplication cp system/fvSolution.SHM system/fvSolution
runApplication decomposePar
runParallel renumberMesh -overwrite
runApplication cp system/controlDict.SRS.init system/controlDict
runApplication cp system/fvSolution.SRS system/fvSolution
runApplication cp system/fvSchemes.SRS system/fvSchemes
runParallel applyBoundaryLayer -ybl 0.1
runApplication cp system/controlDict.SRS.init system/controlDict
runApplication cp system/fvSolution.SRS system/fvSolution
runApplication cp system/fvSchemes.SRS system/fvSchemes
runParallel rhoPimpleFoam
runApplication cp system/controlDict.SRS.avg system/controlDict
runApplication cp system/fvSolution.SRS system/fvSolution
runApplication cp system/fvSchemes.SRS system/fvSchemes
runParallel rhoPimpleFoam
runParallel postProcess -func sampleDict.surface.SRS -latestTime
```

### Step 3: Running the Simulation

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
Now, let's read the `commands.txt` created in [Step 2](#Step-2:-Create-Commands-File).
   ```python
   with open(os.path.join(input_dir,'commands.txt'), 'r') as file:
       commands = [line.strip() for line in file]
   ```
We now have `commands` with a list of commands where each element of that list
is a command.

2. **Run Simulation**:
We now have all we need to run our simulation.
   ```python
   #Choose your simulator
   openfoam = inductiva.simulators.OpenFOAM(distribution="esi")

   task = openfoam.run(
               input_dir=input_dir,
               commands=commands,
               n_vcpus=180,
               use_hwthread=True,
               on=machine_group)
   ```

In this snippet, two arguments might need clarification:

- `n_vcpus`: This sets the number of virtual CPUs (vCPUs) for your simulation,
essentially determining how many parts your simulation will be split into to
run in parallel. Here, we’re dividing the simulation into 180 parts and running
each part simultaneously.

- `use_hwthread`: This enables hyperthreading. Setting this to `True` allows
your simulation to use up to 360 vCPUs on the machine, even if we’re not
utilizing all of them."

3. **Wait and Download Outputs**:
   ```python
   task.wait()
   task.download_outputs()
   ```

4. **Terminate Machine**:
   ```python
   machine_group.terminate()
   ```

### Step 4: Enhancing Performance with MPI Cluster

To further reduce runtime, use an MPI cluster with two machines:

```python
mpi_cluster = inductiva.resources.MPICluster(
                  machine_type="c3d-highcpu-360",
                  data_disk_gb=300,
                  num_machines=2)
mpi_cluster.start()

# Re-run the simulation with adjusted `n_vcpus`
task = openfoam.run(
                  input_dir=input_dir,
                  commands=commands,
                  n_vcpus=360,
                  use_hwthread=True,
                  on=mpi_cluster)
```

**Note**: `spot` machines are not supported for MPI clusters.

### Conclusion

Running the simulation on a high-performance machine and scaling it on an MPI
cluster can significantly reduce computation time. Download and analyze your
results locally once complete. Happy simulating!

## What to read next

Try our [Inductiva API](https://console.inductiva.ai/) 
to streamline your workflows and make the most of cloud resources for 
large-scale simulations.

You may also be interested in reading our blog post,
[The 3D Mesh Resolution Threshold - 5k Points is All You Need!](https://inductiva.ai/blog/article/5k-points-is-all-you-need), 
where we explore just how much you can reduce the level of detail in a 
3D object while still maintaining accurate aerodynamic results in a virtual 
wind tunnel built with OpenFOAM.
