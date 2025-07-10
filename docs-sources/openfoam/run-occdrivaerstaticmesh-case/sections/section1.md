# Prerequisites
Before running the OpenFOAM simulation, make sure all required files are properly set up. This section will walk you through the 
preparation steps.

Let’s get started!

## Download the Simulation Files
Download the case [here](https://develop.openfoam.com/committees/hpc/-/tree/9e0480e778e0c5168b97b8177cc3ece3fb3dc496/incompressible/simpleFoam/occDrivAerStaticMesh) and place it in your working directory under a
folder named `openfoam-occDrivAerStaticMesh/`.

Your directory structure should look like this:

```
- openfoam-occDrivAerStaticMesh/  
  ├── 0.orig
  ├── Allrun
  ├── README.md
  ├── constant
  ├── figures
  └── system
```

The folder structure above represents a typical OpenFOAM case setup:
- `0.orig/` - Contains the initial and boundary conditions for the simulation variables (e.g., velocity, pressure). These files define the starting state of the flow field.
- `constant/` — Stores physical properties of the system such as fluid properties, turbulence model settings, and the mesh definition.
- `system/` — Includes control files that govern the simulation setup, such as solver settings, numerical schemes, and solver algorithms.

This structure is fundamental for OpenFOAM simulations, enabling users to define the problem, specify physical properties, and control numerical execution.

## Download the Mesh Files
Download the 65 million-cell mesh [here](https://zenodo.org/records/15012221/files/polyMesh_65M.tar.gz?download=1).

After downloading, extract the contents into the `constant/` directory. This will create a `polyMesh` folder containing 
all necessary mesh files for the simulation.

## Adapt the `Allrun` script
The provided `Allrun` script is configured for SLURM-based systems, which isn’t applicable when using Inductiva. 
To adapt it, modify the following lines in the `Allrun` file:

```diff
# extra mpi flags
-OMPI_FLAGS="" #  "--mca pml ucx"
+OMPI_FLAGS="--use-hwthread-cpus" #  "--mca pml ucx"

-parEx="mpirun -np ${SLURM_NTASKS} ${OMPI_FLAGS}"
+parEx="mpirun -np ${nProcs} ${OMPI_FLAGS}"
```

## Adjust the Simulation Parameters
By default, this simulation is configured to run with 512 processes, which is quite high for an initial test run. Let’s start with a lower number.

To adjust it, open the file `system/include/caseDefinition` and update the following lines:

```diff
-remover writeInterval_      100;
-nCores              512;              // Number of cores used for simulation
+nCores              180;              // Number of cores used for simulation
decompositionMethod hierarchical;    // Decomposition method
-nHierarchical       (16 8 4);         // Coefficient n for the hierarchical decomposition method
+nHierarchical       (9 5 4);         // Coefficient n for the hierarchical decomposition method
```

This should allow us to run this simulation on a 180 vcpu machine without changing
the grid ratio too much.

That’s it — you’re ready to send your simulation to the cloud!