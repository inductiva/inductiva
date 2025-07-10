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
- Prepare it for execution on an MPI cluster
- Scale the simulation across multiple machines

This guide is divided into a few straightforward steps. Follow along to explore the performance and scalability of this simulation
using the **Inductiva API**.

---

**Prefer to read about it first?** Check out our [blog post](https://inductiva.ai/blog/article/from-supercomputer-to-cloud-a-new-era-for-openfoam-simulations) for the full story.

<div class="cta-bar">
  <div class="cta-text">
    <strong>Ready to dive in?</strong> Click the button to get started with $5 of free credits. No credit card needed!
  </div>
  <button  onclick="window.open('https://console.inductiva.ai/', '_blank')" target="_blank" class="cta-button">Sign In</button>
</div>

```{toctree}
:hidden:
sections/section1.md
sections/section2.md
sections/section3.md
sections/section4.md
```
