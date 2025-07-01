# Prerequisites
Before running the OpenFOAM simulation, ensure all the necessary files are
correctly set up. This guide will walk you through the preparation process.

Let’s get started!

## Download the Simulation File
Download the motorBike example [here](https://develop.openfoam.com/committees/hpc/-/tree/9e0480e778e0c5168b97b8177cc3ece3fb3dc496/incompressible/simpleFoam/occDrivAerStaticMesh) and place it in your working directory under a
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

### Download the Mesh Files

Download the mesh file from [this link](https://zenodo.org/records/15012221/files/polyMesh_65M.tar.gz?download=1)
to obtain the 65 million-cell mesh.
Once downloaded, extract the contents into the `constant/` directory. This will
create a `polyMesh` folder containing all the required mesh files for the simulation.

## Adjust the simulation parameters

By default, this simulation is set to run with 512 processes, which might be a
bit too much for now. To lower this number, open the file
`system/include/caseDefinition` and update the following lines:
```diff
- nCores              512;              // Number of cores used for simulation
+ nCores              112;              // Number of cores used for simulation
decompositionMethod hierarchical;    // Decomposition method
- nHierarchical       (16 8 4);         // Coefficient n for the hierarchical decomposition method
+ nHierarchical       (7 4 4);         // Coefficient n for the hierarchical decomposition method
```

This should allow us to run this simulation on a 112 vcpu machine without changing
the grid ratio too much.

That’s it — you’re ready to send your simulation to the cloud!