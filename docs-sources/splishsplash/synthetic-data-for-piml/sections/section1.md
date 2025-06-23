# Setting Up the Base Case
In the introduction, we outlined the overall workflow for generating synthetic datasets to train Physics-informed ML models using the Inductiva API. Now, let’s dive into the very first step: defining the **base case** simulation model of the physical system we want to study.

Our base case is simple: a **0.5-meter cube of water** is initially positioned in the top corner of a sealed **1-meter cubic box**. At the start of the simulation, the water block is released, causing it to fall, spill, and splash against the closed walls of the box over a 4-second interval. To model this scenario, we use the **SPlisHSPlasH** simulator.

## Preparing the Configuration Files
To get started, we’ve prepared a directory containing all the configuration files needed to run the SPlisHSPlasH simulation. 
Key hyperparameters, most notably the particle radius, have been set to values that enable relatively short simulation times, even 
when using the default computational resources provided by the API.

Let’s begin by [downloading the pre-configured input folder]() and saving it to a local directory. Inside this folder, you’ll find:

- An `.obj` file containing the 3D geometry of the fluid container — in this case, a simple cubic box.  
- A `.json` file specifying the simulation parameters. This configuration file is organized into four key sections that define the simulation:
  - `Configuration`  
  - `RigidBodies`  
  - `Materials`  
  - `FluidModels`

Understanding the structure and contents of the configuration file is essential, as we’ll later modify it programmatically to generate multiple variations of our base case. This will allow us to create a diverse dataset suitable for training a Machine Learning model.

## Configuration File Overview
It's important to understand that these parameters not only govern the physics of the simulation, but also have a significant impact on computational performance. 

```json
{   
    // 1. "Configuration" defines the simulation parameters. For now, we've 
    // modified the particle radius to reduce the simulation runtime.
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
    // 2. "RigidBodies" defines the unit containing the fluid. Here, our
    // RigidBodies are shaped into a cube as outlined in our .obj file
    "RigidBodies": [
        {
            "geometryFile": "unit_box.obj",
            "translation": [0, 0, 0],
            "scale": [1, 1, 1],
            "isDynamic": false
        }
    ],
    // 3. "Materials" defines the properties of the fluid used in the
    // simulation, including its density and viscosity, as well as the
    // viscosity modeling method used by the algorithm.
    "Materials": [
        {
            "id": "Fluid",
            "density0": 1000,
            "viscosity": 1e-6,
            "viscosityMethod": 6
        }
    ],
    // 4. "FluidModels" defines the initial state of the fluid in the
    // simulation. Here, we set a 0.5m fluid cube with no initial velocity.
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

## Running the Base Case
Here is the code required to run this SPlisHSPlasH simulation using the Inductiva API:

```python
import inductiva

# Allocate cloud machine on Google Cloud Platform
cloud_machine = inductiva.resources.MachineGroup(
    provider="GCP",
    machine_type="c2-standard-4")

# Set path to the input directory with the SPlisHSPlasH files
input_dir = "splishsplash-base-dir"

# Initialize the Simulator
splishsplash = inductiva.simulators.SplishSplash()

# Run the simulation
task = splishsplash.run(
    input_dir=input_dir,
    sim_config_filename="config.json",
    on=cloud_machine)

# Wait for the simulation to finish and download the results
task.wait()
cloud_machine.terminate()

task.download_outputs()
```

