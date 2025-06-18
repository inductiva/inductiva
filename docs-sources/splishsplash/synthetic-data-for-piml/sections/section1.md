# Set Up the Base Case 
In the [introduction](synthetic-data-for-piml/index) of this tutorial, we outlined a series of steps required for generating synthetic datasets to train Physics-informed ML models using the Inductiva API. In this section, we’ll dive into the **first step** of that process: defining a base case simulation model of the system under study.

## Base Case Overview
The “base case” simulation is simple: a 0.5m cube of water, initially at rest in one corner of a sealed 1m³ box, is released at the simulation onset. The water is then allowed to splash and interact with the walls of the box over a period of 4 seconds. To simulate this scenario, we’ll use **SPlisHSPlasH**, the SPH simulator used by the authors of the reference study.

> **Note:** If this is your first time using the Inductiva API, we highly recommend setting up your environment via our [User Console](#).

---

## Preparing the Configuration Files

To kick things off, we’ve prepared a pre-configured directory containing all the files needed to run this simulation in SPlisHSPlasH. The key hyperparameter we’ve tuned is the **particle radius**, which directly affects simulation resolution and runtime.

After downloading the input folder locally, you'll find:

- An `.obj` file with the 3D geometry of the fluid container (a simple cubic box).
- A `config.json` file containing simulation parameters. This file consists of four main blocks:

### JSON Configuration Breakdown

```jsonc
{   
    // 1. Simulation configuration settings
    "Configuration": {
        "stopAt": 4,
        "timeStepSize": 0.01,
        "particleRadius": 0.01,
        "simulationMethod": 4,
        "boundaryHandlingMethod": 0,
        "kernel": 1,
        "cflMethod": 1,
        "cflFactor": 0.5,
        "cflMinTimeStepSize": 0.0001,
        "cflMaxTimeStepSize": 0.005,
        "gravitation": [0, 0, -9.81],
        "gradKernel": 1,
        "enableVTKExport": true,
        "dataExportFPS": 60,
        "particleAttributes": "velocity;density"
    },

    // 2. Geometry of the container
    "RigidBodies": [
        {
            "geometryFile": "unit_box.obj",
            "translation": [0, 0, 0],
            "scale": [1, 1, 1],
            "isDynamic": false
        }
    ],

    // 3. Fluid material properties
    "Materials": [
        {
            "id": "Fluid",
            "density0": 1000,
            "viscosity": 1e-6,
            "viscosityMethod": 6
        }
    ],

    // 4. Initial fluid configuration
    "FluidModels": [
        {
            "id": "Fluid",
            "particleFile": "unit_box.obj",
            "translation": [0, 0, 0],
            "scale": [0.5, 0.5, 0.5],
            "initialVelocity": [0, 0, 0]
        }
    ]
}
```

> ⚠️ These parameters control both the **physical behavior** and the **computational cost** of the simulation. In particular, the `particleRadius` has a significant impact on performance and memory usage.

---

## Running the “Base Case”

Let’s now run our base case using the Inductiva API with just a few lines of Python:

```python
import inductiva

# Instantiate machine group
cloud_machine = inductiva.resources.MachineGroup(
    provider="GCP",
    machine_type="c2-standard-4")

# Set path to the input directory
input_dir = "splishsplash-base-dir"

# Initialize the SPlisHSPlasH simulator
splishsplash = inductiva.simulators.SplishSplash()

# Run the simulation
task = splishsplash.run(
    input_dir=input_dir,
    sim_config_filename="config.json",
    on=cloud_machine)

# Wait for completion and download results
task.wait()
cloud_machine.terminate()
task.download_outputs()
```

This script will:
- Upload your local input data to the cloud.
- Schedule and run the simulation.
- Download the output once the task completes.

Terminal output should include the task ID, resource allocation, and progress:

```
Task Information:
> ID:                    i13o8djcoq70bsdut1x73zi69
> Method:                splishsplash
> Local input directory: splishsplash-base-dir
> Submitting to the following computational resources:
 >> Default queue with c2-standard-4 machines.
...
Task i13o8djcoq70bsdut1x73zi69 successfully queued.
```

The simulation will complete in approximately **2 minutes**, and the results will be stored in:

```
inductiva-output/i13o8djcoq70bsdut1x73zi69
├── log
│   └── SPH_log.txt
├── stderr.txt
├── stdout.txt
└── vtk
    ├── ParticleData_Fluid_1.vtk
    ├── ParticleData_Fluid_2.vtk
    ...
```

The `.vtk` files are the **seeds of your dataset**, containing the fluid’s state at each timestep.

---

## Estimating Realistic Computational Costs

Running the above base case was relatively fast, as it only involved ~2,000 particles. In contrast, the study by **Sanchez-Gonzalez et al.** involved simulations with **20,000+ particles**.

To align with that setup, let’s increase the number of particles by decreasing the particle radius from `0.01` to `0.008` in the `Configuration` block:

```jsonc
"Configuration": {
    "stopAt": 4,
    "timeStepSize": 0.01,
    "particleRadius": 0.008
}
```

This change improves resolution but also increases computational cost significantly:
- Runtime: ~30 minutes
- Output size: ~300 MB across 601 frames

If you plan to generate **10,000–20,000 variations**, this would result in **thousands of compute hours** on current hardware.

---

## Scaling Up

To speed things up, you’ll need **more powerful machines** — and fortunately, the Inductiva API makes that incredibly easy. You can adjust compute resources on demand and automate batch jobs to efficiently generate large, high-fidelity datasets for training robust Physics-ML models.

---
