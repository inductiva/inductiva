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

## Example - OpenFOAM Foundation Distribution

```{literalinclude} ../../examples/openfoam-foundation/openfoam-foundation.py
:language: python
```

## Example - ESI Distribution

To run the sample simulation above, simply download the
`openfoam-esi-input-example.zip` file, select the correct distribution by
using `inductiva.simulators.OpenFOAM(distribution="esi")`, and replace  
`runApplication surfaceFeatures` with `runApplication surfaceFeatureExtract`.

## Advanced Example: Running MB9 Micro-benchmark by ExaFOAM

In this advanced example, we’ll be running a complex OpenFOAM simulation 
based on a high-lift configuration setup, which is part of the [ExaFOAM benchmarks](https://exafoam.eu/benchmarks/). Specifically, this case corresponds to the **MB9 micro-benchmark**, which 
serves as preparatory work for the **HPC Grand Challenge** test case of 
the **High Lift Common Research Model (CRM-HL)**. The **CRM-HL** is 
a full aircraft configuration with deployed high-lift devices, simulated 
using **wall-modeled LES (WMLES)**.

The **MB9 micro-benchmark** captures the key characteristics of the **Grand Challenge**, 
such as flow physics and simulation approach, but with significantly fewer 
computational resources. It features a **2D, three-element high-lift wing configuration**, 
simulated with the **IDDES model**, which supports **WMLES** for resolving near-wall 
turbulence.

In this example, we’ll run a modified version of this benchmark using two hardware
configurations powered by **Inductiva**. First, we’ll run the simulation
on a large **360 vCPU machine**, one of the largest single-node configurations
available at Inductiva. Then, to see if we can halve the execution time 
of the simulation, we’ll run the same simulation on an **MPI cluster** using 
two of these large machines.

### Prerequisites

Before running the simulation, make sure to download the input files 
from this [repository](https://develop.openfoam.com/committees/hpc/-/tree/develop/compressible/rhoPimpleFoam/LES/highLiftConfiguration). Once downloaded, place them in a directory named `highLiftConfiguration` 
inside your working folder.

Here’s an example of what your directory structure should look like after 
organizing the necessary files:

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

### Modifications to the Simulation Parameters

To speed up the simulation, we made some adjustments to the original
parameters. Specifically, we set the time step (`dt`) to 0.00002, the
initial time to 0.10, and the final time to 0.30. These changes help reduce the
overall simulation runtime while maintaining reasonable accuracy for this high-lift
configuration case.

### Configuring and Running the Simulation

We'll run the simulation on a **virtual machine with 360 vCPUs** using the
`c3d-highcpu-360` machine type. The process includes several preprocessing
steps, such as **mesh generation** (`blockMesh`, `snappyHexMesh`), followed 
by the simulation run.

The commands used to run the simulation were all taken from the `Allrun` script
included in the downloaded files. We made a few adjustments because all 
commands are executed directly from the `highLiftConfiguration` directory, 
so we can’t use `cd` to navigate into subdirectories and run commands from there.
Additionally, every command must start with `runApplication` or `runParallel`,
even for basic commands like `rm` and `mv`.

We’ve consolidated all the commands into a file called `commands.txt`, with each
command placed on a separate line. This file is then read as a list of strings 
for sequential execution. The contents of `commands.txt` are as follows:

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

Below is the Python code to configure and run the OpenFOAM simulation:

```python notest
import inductiva
import os

machine_group = inductiva.resources.MachineGroup(
    machine_type="c3d-highcpu-360",
    spot=True)
machine_group.start()

# Set simulation input directory
input_dir = "/path/to/highLiftConfiguration"

# Read the simulation commands
with open(os.path.join(input_dir,'commands.txt'), 'r') as file:
    commands = [line.strip() for line in file]

# Initialize the Simulator
openfoam = inductiva.simulators.OpenFOAM(distribution="esi")

# Run simulation
task = openfoam.run(
    input_dir=input_dir,
    commands=commands,
    n_vcpus=180,
    use_hwthread=True,
    on=machine_group)

# Wait for the task to finish and downloading the outputs
task.wait()
task.download_outputs()

# Turn off the machine group
machine_group.terminate()

```

### Important Details

1. **ExaFOAM Benchmark**: This simulation replicates the **MB9 micro-benchmark**
from the [ExaFOAM](https://exafoam.eu/benchmarks/) project. The benchmark
captures key flow physics and simulation strategies of a larger
HPC challenge, but with fewer computational resources, making it scalable 
and efficient for simulating high-lift configurations.
   
2. **Machine Configuration**: We are running the simulation on a `c3d-highcpu-360`
machine with **360 vCPUs** to ensure efficient computation. You can adjust 
the machine configuration based on your available resources or project 
requirements.
   
3. **Commands**: The list of `commands` includes key OpenFOAM tasks
like copying necessary configuration files, running mesh generation, and
parallelizing the process using `decomposePar` and `snappyHexMesh`.

4. **Parallelization**: We utilize `n_vcpus=180` for the simulation, allowing
OpenFOAM to run in parallel and significantly speed up computation. Using 
one thread per physical core, as recommended by [Google Cloud](https://cloud.google.com/blog/products/compute/how-to-reduce-mpi-latency-for-hpc-workloads-on-google-cloud), 
typically offers optimal performance. However, this configuration may not 
always be the best option, as there is no strict rule for selecting the 
number of vCPUs.

### Running the Simulation

Running the script above will take some time, as OpenFOAM tasks like
`snappyHexMesh` and `checkMesh` can be computationally intensive, especially
with large meshes. The `task.wait()` command will block until the simulation
finishes, after which the output files will be available.

In this configuration, the simulation took **40 hours and 20 minutes** to complete,
generating **24.76 GB** of output.

### Post-Simulation

After the simulation completes, the output files will be downloaded locally by
calling `task.download_outputs()`. These outputs include the simulation results,
such as mesh quality reports, final results, and log files for each step in the process.

To inspect the results, navigate to the `inductiva_output` folder, where 
all outputs are stored:

```bash
ls -lasgo inductiva_output/<task_id>
total 672
  0 drwxr-xr-x  143     4576 Sep 19 16:06 .
  0 drwxr-xr-x  400    12800 Sep 19 16:05 ..
  0 drwxr-xr-x   11      352 Sep 19 16:06 0.orig
  8 -rw-r--r--    1      626 Sep 19 16:06 Allclean
 16 -rw-r--r--    1     6998 Sep 19 16:06 Allrun
  8 -rw-r--r--    1      991 Sep 19 16:06 COPYING
 48 -rw-r--r--    1    21547 Sep 19 16:06 README.md
  0 -rw-r--r--    1        0 Sep 19 16:06 case.foam
  0 drwxr-xr-x    8      256 Sep 19 16:06 constant
  0 drwxr-xr-x   27      864 Sep 19 16:06 dynamicCode
  0 drwxr-xr-x   13      416 Sep 19 16:06 figures
  0 drwxr-xr-x    4      128 Sep 19 16:06 processor0
  0 drwxr-xr-x    4      128 Sep 19 16:06 processor1
  0 drwxr-xr-x    4      128 Sep 19 16:06 processor10
    ...
  0 drwxr-xr-x    4      128 Sep 19 16:06 processor126
  0 drwxr-xr-x    4      128 Sep 19 16:06 processor127
 32 -rw-r--r--    1    14437 Sep 19 16:06 stderr.txt
536 -rw-r--r--    1   272331 Sep 19 16:06 stdout.txt
  0 drwxr-xr-x   30      960 Sep 19 16:06 system
 24 -rw-r--r--    1    11399 Sep 19 16:06 thumbnail.png
```

You can now perform post-processing on the results, just as you would if you had
run the simulation locally.

### Utilizing an MPI Cluster for Enhanced Simulation Performance

The initial simulation took a significant amount of time to complete. 
To speed things up, we deployed the simulation across two `c3d-highcpu-360`
machines using an MPI cluster configuration.

Using **Inductiva's Python library**, scaling the simulation resources is
simple. You just need to modify the machine group instantiation to an 
MPI cluster setup and adjust the `n_vcpus` parameter accordingly.

Here’s an example of the code implementation:

```python notest
mpi_cluster = inductiva.resources.MPICluster(
    machine_type="c3d-highcpu-360",
    data_disk_gb=300,
    num_machines=2
)

mpi_cluster.start()

...
# Execute the simulation
task = openfoam.run(
    input_dir=input_dir,
    commands=commands,
    n_vcpus=360, # 360 out of 720 (360 x2)
    use_hwthread=True,
    on=mpi_cluster
)
```

With this simple code change, we were able to significantly reduce the 
simulation time to **19 hours and 10 minutes**, generating an output of **25.02 GB**.

PS: Note that an MPI Cluster does not support `spot` machines.

### Conclusion

This setup shows how cloud resources can be leveraged to run computationally 
intensive **OpenFOAM simulations**, specifically using a case from the **ExaFOAM benchmarks**. 
The cloud-based environment allows for parallel execution on high-performance 
machines, significantly reducing our simulation time by half! Once the simulation 
completes, you can download and analyze the results on your local machine.

Good luck with your OpenFOAM simulations!

## What to read next

If you are interested in OpenFOAM, you may also be interested in checking
the following related simulators that are also available via Inductiva API:

* [CaNS](CaNS.md)
* [DualSPHysics](DualSPHysics.md)
* [SPlisHSPlasH](SPlisHSPlasH.md)

Ready to optimize your own simulations? Try our [Inductiva API](https://console.inductiva.ai/) 
to streamline your workflows and make the most of cloud resources for 
large-scale simulations.

You may also be interested in reading our blog post,
[The 3D Mesh Resolution Threshold - 5k Points is All You Need!](https://inductiva.ai/blog/article/5k-points-is-all-you-need), 
where we explore just how much you can reduce the level of detail in a 
3D object while still maintaining accurate aerodynamic results in a virtual 
wind tunnel built with OpenFOAM.
