
# Scaling it even more

In the previous parts of this tutorial, we progressively scaled our simulation
setup from a single machine to a basic MPI cluster with two nodes. In this part,
we will take the next step by scaling the simulation across **four machines**,
allowing us to further explore how performance scales with available resources.

We'll walk through two different scenarios:

* **With Hyperthreading Enabled:** Where we run the simulation using all available
   virtual CPUs (vCPUs).
* **Without Hyperthreading:** Where we disable hyperthreading to use only the physical cores.

## Scaling to Four Machines

Scaling this simulation to four machines is straightforward. You simply need to
update the `num_machines` parameter in your `MPICluster` configuration to 4 and
you need to edit the `caseDefinition` file to use 448 vCPUs, as shown below:

```diff
-nCores              224;              // Number of cores used for simulation
+nCores              448;              // Number of cores used for simulation
decompositionMethod hierarchical;      // Decomposition method
-nHierarchical       (14 4 4);           // Coefficients for hierarchical decomposition
+nHierarchical       (14 8 4);          // Coefficients for hierarchical decomposition
```

By using 448 vCPUs, this means that we are using hyperthreading, which allows us to
leverage all available virtual CPUs on the four machines.

## Running the Simulation

Here is the code required to run the simulation using the Inductiva API:

```python
import inductiva

# Allocate cloud machine on Google Cloud Platform
cloud_machine = inductiva.resources.MPICluster( \
    provider="GCP",
    machine_type="c2d-highcpu-112",
    data_disk_gb=100,
    num_machines=4,
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

## Disabling Hyperthreading

Turning off hyperthreading is as simple as changing a couple of parameters when
allocatiing your `MPICluster`. You can set the `use_hwthread_cpus` parameter to
`False` and the `np` parameter to `224`, which is the number of physical
cores available across the four machines.

```python
cloud_machine = inductiva.resources.MPICluster( \
    provider="GCP",
    machine_type="c2d-highcpu-112",
    data_disk_gb=100,
    num_machines=4,
    np=224, # Number of processes mpi will use
    use_hwthread_cpus=False, # Disable hyperthreading
    spot=True)
```

Don't forget to revert the `caseDefinition` file to use 224 vCPUs, as shown below:

```diff
nCores              224;              // Number of cores used for simulation
decompositionMethod hierarchical;      // Decomposition method
nHierarchical       (14 4 4);           // Coefficients for hierarchical decomposition
```

If now you run the simulation you will be using only the physical cores available
across the four machines, which is 224 vCPUs.

## Going over the results