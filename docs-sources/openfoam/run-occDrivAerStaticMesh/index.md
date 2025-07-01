# Run occDrivAerStaticMesh from OpenFOAM HPC Benchmark Suite

[This OpenFOAM simulation](https://develop.openfoam.com/committees/hpc/-/tree/9e0480e778e0c5168b97b8177cc3ece3fb3dc496/incompressible/simpleFoam/occDrivAerStaticMesh)
focuses on the open-closed cooling (occ) variant of the
DrivAer model, an industrially relevant benchmark case representing a full-scale
passenger vehicle with static wheels and sealed cooling inlets. The geometry is
based on the notchback version of the Ford Open Cooling DrivAer (OCDA), includes
a highly detailed underbody and engine bay configuration, making it significantly
more complex than the original DrivAer design. For the 2025 OpenFOAM HPC
Challenge, the model was adapted for steady-state RANS simulations using the
SIMPLE algorithm with fixed inner iterations. To ensure compatibility and
scalability in high-performance computing environments, pre-generated meshes of
various resolutions are provided, facilitating robust and reproducible simulations.

This tutorial is split into a few easy steps. Come along and follow through the process with us!


```{toctree}
:hidden:
sections/section1.md
sections/section2.md
sections/section3.md
sections/section4.md
sections/section5.md
```
