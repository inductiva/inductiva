# Run the Simulation on an MPI Cluster
You’re almost ready to scale your simulation across multiple machines, unlocking the full power of an **MPI Cluster**. 
Running simulations in parallel across a cluster can dramatically reduce computation time by distributing the workload over 
many machines.

## Update Simulation Parameters to Use 224 vCPUs
Currently, the simulation is configured to run on 112 vCPUs. To fully utilize the MPI Cluster — which provides 224 vCPUs across two machines — you’ll need to adjust the simulation parameters accordingly.

Open the file `system/include/caseDefinition` and update the number of cores and decomposition settings:

```diff
-nCores              112;              // Number of cores used for simulation
+nCores              224;              // Number of cores used for simulation
decompositionMethod hierarchical;      // Decomposition method
-nHierarchical       (7 4 4);           // Coefficients for hierarchical decomposition
+nHierarchical       (14 4 4);          // Coefficients for hierarchical decomposition
```

These changes ensure that the simulation workload is properly partitioned across all available processors for optimal parallel performance.

## Run your simulation
Below is the updated Python script for running the simulation across two machines using an MPI cluster.

```python
import inductiva

# Allocate cloud machine on Google Cloud Platform
cloud_machine = inductiva.resources.MPICluster( \
    provider="GCP",
    machine_type="c2d-highcpu-112",
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

### How Parallel Commands Are Handled
As mentioned previously, the Inductiva API automatically detects parallel commands by checking for the `-parallel` flag in your command list. When present, it configures the command to run across the entire MPI cluster without requiring any manual setup.

Just like in the single-machine simulation, once the run is complete, the machines are automatically shut down, the results are downloaded, and a detailed summary is printed, as shown below.

```
```

As shown in the “In Progress” line, all the commands (including those executed in parallel) are listed. You can also notice 
that these commands ran on 224 vCPUs, representing the full capacity of our MPI cluster.

With this setup, the computation time was reduced to ---, achieving a --- -fold speedup compared to the previous single-machine simulation.

Keep reading for a deeper dive into performance analysis and various configuration options for your simulation and MPI cluster.



