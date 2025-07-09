# Running Your Simulation on an MPI Cluster

You’re almost ready to scale your simulation across multiple machines, unlocking
the full power of your **MPI Cluster**. Running simulations in parallel on a
cluster allows you to significantly reduce computation time by distributing the
workload across multiple machines.

## Update Simulation Parameters to Use 224 vCPUs

Currently, the simulation is configured to run on 180 vCPUs. To fully utilize
your **MPI Cluster**, which provides 360 vCPUs across multiple machines, you
need to update your configuration accordingly.

Edit the file `system/include/caseDefinition` to increase the number of cores
and adjust the decomposition settings:

```diff
-nCores              180;              // Number of cores used for simulation
+nCores              360;              // Number of cores used for simulation
decompositionMethod hierarchical;      // Decomposition method
-nHierarchical       (9 5 4);           // Coefficients for hierarchical decomposition
+nHierarchical       (15 6 4);          // Coefficients for hierarchical decomposition
```

These changes ensure your simulation workload is properly partitioned across all available processors for optimal parallel performance.

## Run your simulation

Here is the code required to run the simulation using the Inductiva API:

```python
import inductiva

# Allocate cloud machine on Google Cloud Platform
cloud_machine = inductiva.resources.MPICluster( \
    provider="GCP",
    machine_type="c3d-highcpu-180",
    num_machines=2,
    data_disk_gb=100,
    spot=True)

# Initialize OpenFOAM stack
openfoam = inductiva.simulators.OpenFOAM(
    version="2412",
    distribution="esi"
    )

simulation_commands = [
    "cp system/controlDict.noWrite system/controlDict",
    "cp system/fvSolution.fixedIter system/fvSolution",
    "decomposePar -constant",
    "restore0Dir -processor",
    "renumberMesh -constant -overwrite -parallel",
    "potentialFoam -initialiseUBCs -parallel",
    "applyBoundaryLayer -ybl '0.0450244' -parallel",
    "simpleFoam -parallel",
]

task = openfoam.run( \
    input_dir="/Path/to/openfoam-occDrivAerStaticMesh",
    commands=simulation_commands,
    on=cloud_machine)

# Wait for the simulation to finish and download the results
task.wait()
cloud_machine.terminate()

task.download_outputs()

task.print_summary()
```

### Explaining the parallel commands

Like we said before, we automatically detect parallel commands and configure them
to run on your MPI Cluster. We look for the flag `-parallel` in the command
list, and if it is present, we automatically do the necessary configuration.

Just like the previous simulation, when the simulation is complete, we terminate
the machine, download the results and print a summary of the simulation as
shown.

```
```

As you can see in the “In Progress” line, the part of the timeline that
represents the actual execution of the simulation, we have all our commands that
ran during the simulation, including the parallel commands. You can also see
that those commands ran on 360 vCPUs, which is the total number of vCPUs available
on our MPI Cluster.

With this configuration, we managed to reduce the computaion time to X time, meaning
a reduction of X times compared to the previous single machine simulation.

Keep reading for a more in depth dive into performance and different configurations
of your simulation and MPI Cluster.