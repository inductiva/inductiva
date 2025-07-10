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
  <button class="cta-button">Sign In</button>
</div>

.cta-bar {
  display: flex;
  justify-content: space-between;
  align-items: center;
  gap: 16px;
  background-color: #f3f4fe;
  padding: 20px 28px;
  border-radius: 10px;
  box-shadow: 0 4px 12px rgba(0, 0, 0, 0.06);
  max-width: 800px;
  margin: 40px auto;
  flex-wrap: wrap; /* responsive */
}

.cta-text {
  font-size: 16px;
  color: #333;
  flex: 1;
  min-width: 250px;
}

.cta-button {
  background-color: #8E3BFF; /* your original purple */
  color: white;
  border: none;
  padding: 12px 24px;
  border-radius: 8px;
  font-size: 16px;
  font-weight: bold;
  cursor: pointer;
  transition: background-color 0.3s ease;
}

.cta-button:hover {
  background-color: #752fd6; /* darker shade of your purple */
}