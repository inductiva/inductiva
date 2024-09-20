# OpenFOAM

OpenFOAM is a Finite Volume method for CFD simulations with a wide range of 
applications across several areas of engineering and science. OpenFOAM has an 
extensive
range of features to solve anything from complex fluid flows involving chemical 
reactions, turbulence and heat transfer, to solid dynamics and electromagnetics.

There are two main open-source distributions of OpenFOAM, one developed by
[OpenFOAM foundation](https://openfoam.org/) and another by the
[ESI Group](https://www.openfoam.com/). Inductiva API supports both distributions,
and you can choose which one you want to use by setting the `distribution` parameter
when initializing the simulator. The default distribution is the one developed by
the OpenFOAM foundation.

A single simulation via Inductiva API comprises several steps done via OpenFOAM 
- e.g., partitioning the domain, meshing, solvers and post-processing. Hence, to 
configure a simulation for OpenFOAM the user will need a set of configuration
files that are organized in three folders:
- `time`: containing individual files of data for particular fields, like
initial values and boundary conditions that the user must specify. For example,
for an initial condition at $, the initial conditions will be stored in the
directory `0`.
- `constant`: contains files that describe the objects in the simulation and the 
physical properties of the application we are concerned.
- `system`: contains all of the files that describe the simulation, including
the solvers, the numerical parameters, and the output files. It must contain at
least 3 files: `controlDict` where run control parameters are set including
start/end time, time step and parameters for data output; `fvSchemes` where
discretization schemes used in the solution may be selected at run-time; and
`fvSolution` where the equation solvers, tolerances and other algorithm controls
are set for the run.

All of these folders should be inside an input directory. Finally, to run a 
simulation the user needs to configure a list of dictionaries specifying the
commands they want to execute on the backend. Below, we run the
[motorbike tutorial](https://github.com/OpenFOAM/OpenFOAM-8/tree/master/tutorials/incompressible/simpleFoam/motorBike) from OpenFOAM and show how this is done in practice.

The commands passed to the simulator follow the structure of OpenFOAM, that is, 
using the prefix `runApplication` the command will execute sequentially and with
`runParallel` the command will use the maximum number of cores the machine has 
available. Hence, you don't need to set the specific number of processes, the 
simulator will do that for you automatically. In particular, the
decomposeParDict will be configured automatically and, at the moment, only the 
scotch decomposition method is available.

## Example - Foundation distribution

````python
import inductiva

# Instantiate machine group
machine_group = inductiva.resources.MachineGroup('c2-standard-4')
machine_group.start()

# Set simulation input directory
input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "openfoam-input-example.zip", unzip=True)

# Set the simulation commands
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

# Initialize the Simulator
openfoam = inductiva.simulators.OpenFOAM(distribution="foundation")

# Run simulation with config files in the input directory
task = openfoam.run(input_dir=input_dir,
                    commands=commands,
                    n_vcpus=4,
                    on=machine_group)

task.wait()
task.download_outputs()

machine_group.terminate()
````

## Example - ESI distribution

````python
import inductiva

# Instantiate machine group
machine_group = inductiva.resources.MachineGroup('c2-standard-4')
machine_group.start()

# Set simulation input directory
input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "openfoam-esi-input-example.zip", unzip=True)

# Set the simulation commands
commands = [
    "runApplication surfaceFeatureExtract",
    "runApplication blockMesh",
    "runApplication decomposePar -copyZero",
    "runParallel snappyHexMesh -overwrite",
    "runParallel potentialFoam",
    "runParallel simpleFoam",
    "runApplication reconstructParMesh -constant",
    "runApplication reconstructPar -latestTime"
]

# Initialize the Simulator
openfoam = inductiva.simulators.OpenFOAM(distribution="esi")

# Run simulation with config files in the input directory
task = openfoam.run(input_dir=input_dir,
                    commands=commands,
                    n_vcpus=4,
                    on=machine_group)

task.wait()
task.download_outputs()

machine_group.terminate()
````

## Advanced Example: OpenFOAM Simulation (ExaFOAM Benchmark)

Here's a guide for running an OpenFOAM simulation using the Inductiva platform,
focusing on a high-lift configuration case from the `rhoPimpleFoam` solver.
This specific example is based on one of the benchmark cases provided by the
[ExaFOAM project](https://exafoam.eu/benchmarks/), which aims to develop highly
scalable CFD solvers using OpenFOAM for exascale computing. We'll use a powerful
cloud setup to handle the computationally intensive mesh generation and simulation.

We will run a complex OpenFOAM simulation based on a high-lift configuration
setup, which is part of the [ExaFOAM benchmarks](https://exafoam.eu/benchmarks/).
Specifically, this case corresponds to the MB9 micro-benchmark, which is
preparatory work leading up to the HPC Grand Challenge test case of the High
Lift Common Research Model (CRM-HL). The CRM-HL is a full aircraft configuration
with deployed high-lift devices, simulated using wall-modeled LES (WMLES). 

The MB9 micro-benchmark captures the essential characteristics of the Grand
Challenge (such as flow physics and simulation approach) while requiring
significantly fewer computational resources. It features a two-dimensional,
three-element high-lift wing configuration, simulated with the IDDES model,
which provides WMLES functionality in regions of resolved near-wall turbulence.
This case is based on the well-known 30P30N test case, extensively studied in
the 4th AIAA CFD High Lift Prediction Workshop (HLPW-4) and supported by
available experimental data.

You can download the input files for this specific case from the
[OpenFOAM HPC repository](https://develop.openfoam.com/committees/hpc/-/tree/develop/compressible/rhoPimpleFoam/LES/highLiftConfiguration).

### Prerequisites

Before running the simulation, ensure you have downloaded the input files. You
can get them from the following
[repository](https://develop.openfoam.com/committees/hpc/-/tree/develop/compressible/rhoPimpleFoam/LES/highLiftConfiguration).
Once downloaded, move them into a directory named `highLiftConfiguration` inside
your working folder.

Here’s an example of how your directory structure should look after downloading
the necessary files:

```bash
ls -lasgo highLiftConfiguration
total 104
 0 drwxrwxr-x@  12     384 Jun 13 10:09 .
 0 drwx------@ 234    7488 Sep 19 15:15 ..
 0 drwxrwxr-x@  11     352 Jun 13 10:09 0.orig
 8 -rwxr-xr-x@   1     626 Jun 13 10:09 Allclean
16 -rwxr-xr-x@   1    6998 Jun 13 10:09 Allrun
 8 -rw-rw-r--@   1     991 Jun 13 10:09 COPYING
48 -rw-rw-r--@   1   21547 Jun 13 10:09 README.md
 0 -rw-rw-r--@   1       0 Jun 13 10:09 case.foam
 0 drwxrwxr-x@   5     160 Jun 13 10:09 constant
 0 drwxrwxr-x@  13     416 Jun 13 10:09 figures
 0 drwxrwxr-x@  27     864 Jun 13 10:09 system
24 -rw-rw-r--@   1   11399 Jun 13 10:09 thumbnail.png
```

### Configuring and Running the Simulation

We'll run the simulation on a virtual machine with 180 CPUs using the
`c3d-standard-180` machine type. The simulation involves several preprocessing
steps, including mesh generation (`blockMesh`, `snappyHexMesh`), and finally,
the simulation run itself. 

Below is the Python code to configure and run the OpenFOAM simulation:

```python
import inductiva
import os

machine_group = inductiva.resources.MachineGroup(machine_type="c3d-standard-180")
machine_group.start()

import inductiva

# Set simulation input directory
input_dir = "/path/to/highLiftConfiguration"

# Read the simulation commands
with open(os.path.join(input_dir,'input.txt'), 'r') as file:
    commands = [line.strip() for line in file]

# Initialize the Simulator
openfoam = inductiva.simulators.OpenFOAM(distribution="esi")

# Run simulation
task = openfoam.run(
    input_dir=input_dir,
    commands=commands,
    n_vcpus=128,
    use_hwthread=True,
    on=machine_group)

# Wait for the task to finish and downloading the outputs
task.wait()
task.download_outputs()

# Turn off the machine group
machine_group.terminate()

```

The commands used in this simulation were taken from the `Allrun` script found
in the downloaded files. Some adjustments were made because we execute all
commands directly from the `highLiftConfiguration` directory, meaning we can't
use `cd` to navigate into subdirectories and run commands from there.
Additionally, every command must start with `runApplication` or `runParallel`,
which is why even basic commands like `rm` and `mv` are preceded by `runApplication`.

We have consolidated all the commands into a file named `input.txt`, with each
command placed on a separate line. This file is then read into a list of strings
to execute sequentially. The contents of `input.txt` are as follows:
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
runParallel applyBoundaryLayer -ybl 0.1 > ${LOGDIR}/log.R03.applyBoundaryLayer 2>&1 || exit 1
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

### Important Details

1. **ExaFOAM Benchmark**: This simulation replicates the MB9 micro-benchmark
from the [ExaFOAM](https://exafoam.eu/benchmarks/) project. The benchmark
focuses on capturing key flow physics and the simulation approach of a larger
HPC challenge with fewer computational resources, providing a scalable and
efficient simulation of a high-lift configuration.
   
2. **Machine Configuration**: We are using a powerful `c3d-standard-180`
machine with 180 CPUs to ensure that the computation proceeds efficiently.
You can adjust this depending on the available resources or your project's needs.
   
3. **Commands**: The list of `commands` involves several key OpenFOAM tasks
like copying necessary configuration files, running mesh generation, and
parallelizing the process using `decomposePar` and `snappyHexMesh`.

4. **Parallelization**: We use `n_vcpus=128` in the simulation, which will
allow OpenFOAM to run in parallel, greatly speeding up the computation.

### Running the Simulation

Running the script above will take some time, as OpenFOAM tasks like
`snappyHexMesh` and `checkMesh` can be computationally intensive, especially
with large meshes. The `task.wait()` command will block until the simulation
is complete, after which the output files will be available.

### Post-Simulation

After the simulation completes, the output files will be downloaded locally by
calling `task.download_outputs()`. The outputs will include simulation results
like mesh quality reports, final results, and log files for each step in the process.

To inspect the output, look under the `inductiva_output` folder where the
results are downloaded:

```bash
ls -lasgo inductiva_output/task_id
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

### Small benchmark

We conducted a series of OpenFOAM simulations across different machine
configurations to evaluate the performance impact of varying core counts and the
effect of hyperthreading. The configurations ranged from 90 to 360 virtual CPUs
(vCPUs), with some tests utilizing hyperthreading while others disabled it. Our
goal was to assess how these hardware differences influence simulation time.

| Machine Configuration            | nCores/n_vCPUs           | Hyperthreading Status        | Time (min:sec) |
|-----------------------------------|--------------------------|------------------------------|----------------|
| c3d-standard-90                   | 90                       | Enabled                      | 7:48           |
| c3d-standard-180                  | 180                      | Enabled                      | 7:39           |
| c3d-standard-360                  | 360                      | Enabled                      | 8:32           |
| c3d-standard-90                   | 45                       | Disabled (Hyperthreading Off) | 7:20           |
| c3d-standard-180                  | 90                       | Disabled (Hyperthreading Off) | 6:46           |
| c3d-standard-360                  | 180                      | Disabled (Hyperthreading Off) | 6:57           |


The results clearly demonstrate that disabling hyperthreading significantly
enhances simulation performance across all tested configurations. For instance,
in the c3d-standard-180 setup, disabling hyperthreading reduced the runtime by
nearly a minute, from 7:39 to 6:46. However, the c3d-standard-360 configuration,
with 360 cores, performed slower than the c3d-standard-180 machines, taking 8:32
to complete. This suggests that the current OpenFOAM simulation may not scale
efficiently across such a large number of cores, likely due to the overhead of
managing many processes relative to the actual work performed by each. While
larger machines like the c3d-standard-360 may not be ideal for this particular
simulation, they could still offer advantages for other types of simulations
better suited to higher core counts.

### Conclusion

This setup demonstrates how to leverage cloud resources to run computationally
heavy OpenFOAM simulations, specifically using a case from the ExaFOAM benchmarks.
The cloud-based environment enables parallel execution on high-performance
machines, drastically reducing computation time. Once the simulation completes,
you can download and analyze the results on your local machine.

Good luck with your OpenFOAM simulations!

## What to read next

If you are interested in OpenFOAM, you may also be interested in checking
the following related simulators that are also available via Inductiva API:

* [CaNS](CaNS.md)
* [DualSPHysics](DualSPHysics.md)
* [SPlisHSPlasH](SPlisHSPlasH.md)

You may also be interested in reading our blog post
[The 3D Mesh Resolution Threshold - 5k Points is All You Need!](https://inductiva.ai/blog/article/5k-points-is-all-you-need),
where we investigate the impact of reducing the level of detail of a 3D object in
the accuracy of aerodynamic metrics obtained using a (virtual) wind tunnel
implemented using OpenFOAM.
