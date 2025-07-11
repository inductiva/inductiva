# Prerequisites
Before running an OpenFOAM simulation, ensure all the necessary files are correctly set up. This guide will walk you through the preparation process.

Let’s get started!

<p align="center"><img src="../../_static/bike_flow_slice_U.png" alt="OpenFOAM simulation visualization" width="700"></p>

## Download the Simulation Files
Download the motorBike example [here](https://storage.googleapis.com/inductiva-api-demo-files/openfoam-input-example.zip) and place it in your working directory under a folder named `openfoam-input-example/`.

Your directory structure should look like this:

```
- openfoam-input-example/  
  ├── 0
  ├── Allrun
  ├── constant
  └── system
```

The folder structure above represents a typical OpenFOAM case setup:
- `0/` - Contains the initial and boundary conditions for the simulation variables (e.g., velocity, pressure). These files define the starting state of the flow field.
- `constant/` — Stores physical properties of the system such as fluid properties, turbulence model settings, and the mesh definition.
- `system/` — Includes control files that govern the simulation setup, such as solver settings, numerical schemes, and solver algorithms.

This structure is fundamental for OpenFOAM simulations, enabling users to define the problem, specify physical properties, and control numerical execution.

That’s it — you’re ready to send your simulation to the cloud!