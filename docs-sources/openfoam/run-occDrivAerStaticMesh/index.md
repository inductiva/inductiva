# Run `occDrivAerStaticMesh` from the OpenFOAM HPC Benchmark Suite
The [occDrivAerStaticMesh simulation](https://develop.openfoam.com/committees/hpc/-/tree/9e0480e778e0c5168b97b8177cc3ece3fb3dc496/incompressible/simpleFoam/occDrivAerStaticMesh) focuses on the open-closed cooling (occ) variant of the DrivAer model - 
an industrially relevant benchmark case that represents a full-scale passenger vehicle with static wheels and sealed cooling inlets.

The geometry is derived from the notchback version of the **Ford Open Cooling DrivAer (OCDA) model** and features a highly 
detailed underbody and engine bay configuration, adding significant complexity compared to the original DrivAer design.

As part of the [2025 OpenFOAM HPC Challenge](https://wiki.openfoam.com/images/4/44/OHC-1.pdf), the model was adapted for steady-state 
RANS simulations using the SIMPLE algorithm with fixed inner iterations. To support compatibility and scalability in 
high-performance computing environments, pre-generated meshes at multiple resolutions are provided, enabling robust and 
reproducible simulations.

In this tutorial, we'll walk through how to:
- Run the `occDrivAerStaticMesh` simulation
- Scale the simulation across multiple machines using an **MPI cluster**
- Benchmark performance with hyperthreading enabled and disabled

This guide is divided into a few straightforward steps. Follow along to explore the performance and scalability of this simulation 
using the **Inductiva API**!

```{toctree}
:hidden:
sections/section1.md
sections/section2.md
sections/section3.md
sections/section4.md
sections/section5.md
```
