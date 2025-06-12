# Run the 30P30N micro-benchmark (MB9) from ExaFOAM
In this tutorial we will show you how to use the Inductiva API to run an
advanced OpenFOAM case that requires significant computing power.

## Objective
The goal of this tutorial is to demonstrate how to run the
`30P30N micro-benchmark from ExaFOAM` use case from the CFD Tutorials, available in the [official ExaFOAM documentation](https://exafoam.eu/benchmarks/).

## Prerequisites
Download the required files [here](https://develop.openfoam.com/committees/hpc/-/tree/6b9c8438ec54961d2968977f84a47a5f9c5be4b6/compressible/rhoPimpleFoam/LES/highLiftConfiguration) and place them in a folder named `highLiftConfiguration`.

The **directory structure** should look like this:
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
 
## Running Your Simulation
Here is the code required to run a OpenFOAM simulation using the Inductiva API:

```python
"""OpenFOAM example"""
import inductiva

# Allocate cloud machine on Google Cloud Platform
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c3d-standard-360",
	spot=True)

# Initialize the Simulator
OpenFOAM = inductiva.simulators.OpenFOAM( \
    version="2412",
	distribution="esi")

# Run simulation
task = OpenFOAM.run(input_dir="/path/to/highLiftConfiguration",
    shell_script="./Allrun",
    on=cloud_machine)

# Wait for the simulation to finish and download the results
task.wait()
cloud_machine.terminate()

task.download_outputs()

task.print_summary()
```

When the simulation is complete, we terminate the machine, download the results and print a summary of the simulation as shown below.

```
asd
```

As you can see in the "In Progress" line, the part of the timeline that represents the actual execution of the simulation, the core computation time 
of this simulation was approximately .

## Scalling your simulation to multiple machines

As you can see above, some simulations can take a long time to finish. For the more
intensive simulations we offer our users the ability to use MPI clusters. MPI clusters
are a type of resource with multiple machines that behaves as a single very powerfull
machines. This means, if you have an MPI cluster with 2 machines with 100 cores,
you can run your simulation with a total of 200 cores.

In order to take advantage of this feature we can't run our simulation using the
`shell_script` argument. We need to run our simulation with a new argument, `commands`.

This argument will have a list off all commands that our simulation will run. Some of this
commands will run sequencialy, and other commands will run in parallel.

### Result of converting the Allrun into a list of commands

For this part we just need to convert all commands inside the shell script into
a list of single commands. We just have one caveat. Commands that will run in parallel
need to be defined with our MPIConfig and Command class (as you can see below).

```python
#Mpi configuration to use for the parallel commands
mpiconfig = inductiva.commands.MPIConfig("4.1.6",np=1344,use_hwthread_cpus=True)

#Parallel commands for the simulation
command1 = inductiva.commands.Command("snappyHexMesh -parallel -decomposeParDict system/decomposeParDict.SHM -dict system/snappyHexMeshDict -overwrite",mpi_config=mpiconfig)
command2 = inductiva.commands.Command("checkMesh -parallel -decomposeParDict system/decomposeParDict.SHM -latestTime -meshQuality",mpi_config=mpiconfig)
command3 = inductiva.commands.Command("renumberMesh -parallel -overwrite",mpi_config=mpiconfig)
command4 = inductiva.commands.Command("applyBoundaryLayer -parallel -ybl 0.1",mpi_config=mpiconfig)    
command5 = inductiva.commands.Command("rhoPimpleFoam -parallel",mpi_config=mpiconfig)
command6 = inductiva.commands.Command("rhoPimpleFoam -parallel",mpi_config=mpiconfig)
command7 = inductiva.commands.Command("postProcess -func sampleDict.surface.SRS -parallel -latestTime",mpi_config=mpiconfig)

#List of all commands, sequential and parallel
commands = ["cp system/controlDict.SHM system/controlDict",
            "cp system/fvSchemes.SHM system/fvSchemes",
            "cp system/fvSolution.SHM system/fvSolution",
            "blockMesh",
            "snappyHexMesh -dict system/snappyHexMeshDict.refineblockMesh -overwrite",
            "mv ./0/cellLevel ./constant/polyMesh/",
            "mv ./0/pointLevel ./constant/polyMesh/",
            "rm -r 0",
            "checkMesh",
            "extrudeMesh -dict system/extrudeMeshDict.refineblockMesh",
            "rm constant/polyMesh/cellLevel",
            "rm constant/polyMesh/cellZones",
            "rm constant/polyMesh/pointLevel",
            "rm constant/polyMesh/pointZones",
            "rm constant/polyMesh/faceZones",
            "rm constant/polyMesh/level0Edge",
            "rm constant/polyMesh/surfaceIndex",
            'sed -i s/"ff_zMin"/"ff_SAVE"/g constant/polyMesh/boundary',
            'sed -i s/"ff_zMax"/"ff_zMin"/g constant/polyMesh/boundary',
            'sed -i s/"ff_SAVE"/"ff_zMax"/g constant/polyMesh/boundary',
            "checkMesh",
            "surfaceFeatureExtract",
            "decomposePar -force -decomposeParDict system/decomposeParDict.SHM",
            command1,
            command2,
            "reconstructParMesh -mergeTol 1e-08 -constant -latestTime",
            "echo -e -n '\n    Mirror 3D mesh to apply cyclic BC ...\n'",
            "cp -r constant/polyMesh constant/polyMesh.origHalf",
            "mirrorMesh -overwrite",
            "topoSet -dict system/topoSetDict.faces.ff_zMin",
            "createPatch -overwrite -dict system/createPatchDict.mirrorMesh",
            "changeDictionary -constant -dict system/changeDictionaryDict.cyclicPatches -enableFunctionEntries",
            "topoSet -dict system/topoSetDict.faces.cyclic",
            "checkMesh -constant",
            "rm -r constant/polyMesh.origHalf",
            "cp -r 0.orig 0",
            "cp system/controlDict.SHM system/controlDict",
            "cp system/fvSchemes.SHM system/fvSchemes",
            "cp system/fvSolution.SHM system/fvSolution",
            "decomposePar -force",
            command3,
            "cp system/controlDict.SRS.init system/controlDict",
            "cp system/fvSolution.SRS system/fvSolution",
            "cp system/fvSchemes.SRS system/fvSchemes",
            command4,
            "echo -e -n '\n    Running initial SRS ...\n'",
	        "cp system/controlDict.SRS.init system/controlDict",
            "cp system/fvSolution.SRS system/fvSolution",
            "cp system/fvSchemes.SRS system/fvSchemes",
            command5,
            "echo -e -n '\n    Running productive SRS ...\n'",
            "cp system/controlDict.SRS.avg system/controlDict",
            "cp system/fvSolution.SRS system/fvSolution",
            "cp system/fvSchemes.SRS system/fvSchemes",
            command6,
            command7,
            ]
```

Once you have this change applies, running the simulation should now be simple.

## Running Your Simulation on an MPI cluster

```python
#MPI CLUSTER asd
"""OpenFOAM example"""
import inductiva

# Allocate cloud machine on Google Cloud Platform
cloud_machine = inductiva.resources.MPICluster(
    machine_type="c2d-highcpu-112",
    data_disk_gb=300,
    num_machines=12,
    spot=True)


# Initialize the Simulator
OpenFOAM = inductiva.simulators.OpenFOAM( \
    version="2412",
	distribution="esi")

mpiconfig = inductiva.commands.MPIConfig("4.1.6",np=1344,use_hwthread_cpus=True)

command1 = inductiva.commands.Command("snappyHexMesh -parallel -decomposeParDict system/decomposeParDict.SHM -dict system/snappyHexMeshDict -overwrite",mpi_config=mpiconfig)
command2 = inductiva.commands.Command("checkMesh -parallel -decomposeParDict system/decomposeParDict.SHM -latestTime -meshQuality",mpi_config=mpiconfig)
command3 = inductiva.commands.Command("renumberMesh -parallel -overwrite",mpi_config=mpiconfig)
command4 = inductiva.commands.Command("applyBoundaryLayer -parallel -ybl 0.1",mpi_config=mpiconfig)    
command5 = inductiva.commands.Command("rhoPimpleFoam -parallel",mpi_config=mpiconfig)
command6 = inductiva.commands.Command("rhoPimpleFoam -parallel",mpi_config=mpiconfig)
command7 = inductiva.commands.Command("postProcess -func sampleDict.surface.SRS -parallel -latestTime",mpi_config=mpiconfig)

commands = ["cp system/controlDict.SHM system/controlDict",
            "cp system/fvSchemes.SHM system/fvSchemes",
            "cp system/fvSolution.SHM system/fvSolution",
            "blockMesh",
            "snappyHexMesh -dict system/snappyHexMeshDict.refineblockMesh -overwrite",
            "mv ./0/cellLevel ./constant/polyMesh/",
            "mv ./0/pointLevel ./constant/polyMesh/",
            "rm -r 0",
            "checkMesh",
            "extrudeMesh -dict system/extrudeMeshDict.refineblockMesh",
            "rm constant/polyMesh/cellLevel",
            "rm constant/polyMesh/cellZones",
            "rm constant/polyMesh/pointLevel",
            "rm constant/polyMesh/pointZones",
            "rm constant/polyMesh/faceZones",
            "rm constant/polyMesh/level0Edge",
            "rm constant/polyMesh/surfaceIndex",
            'sed -i s/"ff_zMin"/"ff_SAVE"/g constant/polyMesh/boundary',
            'sed -i s/"ff_zMax"/"ff_zMin"/g constant/polyMesh/boundary',
            'sed -i s/"ff_SAVE"/"ff_zMax"/g constant/polyMesh/boundary',
            "checkMesh",
            "surfaceFeatureExtract",
            "decomposePar -force -decomposeParDict system/decomposeParDict.SHM",
            command1,
            command2,
            "reconstructParMesh -mergeTol 1e-08 -constant -latestTime",
            "echo -e -n '\n    Mirror 3D mesh to apply cyclic BC ...\n'",
            "cp -r constant/polyMesh constant/polyMesh.origHalf",
            "mirrorMesh -overwrite",
            "topoSet -dict system/topoSetDict.faces.ff_zMin",
            "createPatch -overwrite -dict system/createPatchDict.mirrorMesh",
            "changeDictionary -constant -dict system/changeDictionaryDict.cyclicPatches -enableFunctionEntries",
            "topoSet -dict system/topoSetDict.faces.cyclic",
            "checkMesh -constant",
            "rm -r constant/polyMesh.origHalf",
            "cp -r 0.orig 0",
            "cp system/controlDict.SHM system/controlDict",
            "cp system/fvSchemes.SHM system/fvSchemes",
            "cp system/fvSolution.SHM system/fvSolution",
            "decomposePar -force",
            command3,
            "cp system/controlDict.SRS.init system/controlDict",
            "cp system/fvSolution.SRS system/fvSolution",
            "cp system/fvSchemes.SRS system/fvSchemes",
            command4,
            "echo -e -n '\n    Running initial SRS ...\n'",
	        "cp system/controlDict.SRS.init system/controlDict",
            "cp system/fvSolution.SRS system/fvSolution",
            "cp system/fvSchemes.SRS system/fvSchemes",
            command5,
            "echo -e -n '\n    Running productive SRS ...\n'",
            "cp system/controlDict.SRS.avg system/controlDict",
            "cp system/fvSolution.SRS system/fvSolution",
            "cp system/fvSchemes.SRS system/fvSchemes",
            command6,
            command7,
            ]

# Run simulation
task = OpenFOAM.run(
    input_dir="/path/to/highLiftConfiguration",
    commands=commands,
    on=cloud_machine)

# Wait for the simulation to finish and download the results
task.wait()
cloud_machine.terminate()

task.download_outputs()

task.print_summary()

```

Here is the result of running this simulation on a couple of MPI cluster configurations:

| Machine Type       | Num Machines | Cores | Duration         | Cost             |
|--------------------|--------------|-------|------------------|------------------|
| c3d-highcpu-360    | 1            | 360   | 41 h             | 140.35 US$       |
| c3d-highcpu-360    | 3            | 1080  | 17h 11min        | 179.23 US$       |
| c2d-highcpu-112    | 8            | 896   | 14h 59min        | 92.29 US$        |
| c2d-highcpu-112    | 12           | 1344  | 12h 5min (Preempted after 12h 5min) | N/A |

As you can see, we can improve the simulation time a lot, but we can't scale forever.
In the last example we went from 896 cores to 1344 and we were almost having the same
duration (at least the same, it could even be slower) as the one with 896 cores.

Feel free to explore with new MPI cluster combinations and happy simulations!