---
myst:
  html_meta:
    description: "Dive into the first crucial step of synthetic data generation and learn how to define your 'base case' simulation model."
    keywords: "Inductiva API, Programming, HPC, Simulation, Tutorial, Synthetic Data Generation, Physics-ML, SPH"
---

# Setting Up the "Base Case"

## Quick Recap

In the [introduction](synthetic-data-generation-1.md) of this tutorial, we outlined a series of steps needed for generating synthetic datasets for training Physics-ML models using the Inductiva API. So, in this chapter, we will dive into the first step of this process: **defining a "base case"** simulation model of the system we wish to study.

The "base case" simulation is quite simple: a 0.5m cube of water, initially at rest at one of the top corners of a sealed 1m cubic box, is dropped at the simulation onset, allowing the water to spill and splash against the walls of the closed box 
for 4 seconds. For simulating this base case, we will be using [SPlisHSPlasH](https://tutorials.inductiva.ai/simulators/SPlisHSPlasH.html), the SPH simulator used by the authors.

> Note: if this is your first encounter with our 
API, we highly recommend going through our [User Console](https://console.inductiva.ai/) to set up your environment correctly. 

## Preparing the Configuration Files

To kick things off, we've pre-configured a directory containing all the configuration files necessary to run the SPlisHSPlasH simulation. We defined relevant hyperparameters, namely the particle radius, with values that allow for relatively short simulation times, even using the default computational resources available via the API. 

>Let's **<a href="../_static/generating-synthetic-data/splishsplash-base-dir.zip" download="splishsplash-base-dir.zip" class="bi bi-cloud-download-fill">download our pre-configured input folder, </a>** and store it in a local directory.
      
In this folder, we'll find:

- An `.obj` file with the 3D geometry of 
the fluid container, in this case a simple  cubic box.

- A `JSON` file containing the simulation parameters. This file essentially comprises **four key blocks** that define the whole simulation; 
`Configuration`, `RigidBodies`, `Materials` and `FluidModels`. Understanding the contents of this cofiguration file is crucial because we will later  
tweak it to be able to programmatically produce multiple variations of our "base case" and generate 
a diverse enough dataset to train a machine learning (ML) model.

Let's take a closer look at these blocks:

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
It's important to understand that these parameters not only control the physics of the simulation, but have also profound impact in how computation is made. Most specficially, and as we will see later, the **particle size** plays a crucial role in determining the computational cost of the simulation.

## Running our "Base Case"

To run our "base case", we'll initialize the SPlisHSPlasH simulator using the API 
with just a couple lines of code:

```python
import inductiva

# Instantiate machine group
cloud_machine = inductiva.resources.MachineGroup(
    provider="GCP",
    machine_type="c2-standard-4")

# Set path to the input directory with the SPlisHSPlasH configuration files
input_dir = "splishsplash-base-dir"

# Initialize the SPlisHSPlasH simulator through the API
splishsplash = inductiva.simulators.SplishSplash()

# Run the simulation task with the parameters defined in the .json file
task = splishsplash.run(input_dir=input_dir,
                        sim_config_filename="config.json",
                        on=cloud_machine)

# Wait for the simulation to complete and download the outputs
task.wait()
# Terminate the machine group
cloud_machine.terminate()

task.download_outputs()
```
This script will upload the input data from our local directory to the API server and schedule a simulation `task` for execution. We will be able check details about the `task`, including its ID and the machine group assigned for its computation, by observing the stdout of your terminal:

```
Task Information:
> ID:                    i13o8djcoq70bsdut1x73zi69
> Method:                splishsplash
> Local input directory: splishsplash-base-dir
> Submitting to the following computational resources:
 >> Default queue with c2-standard-4 machines.
Preparing upload of the local input directory splishsplash-input-example
(128 B).
Local input directory successfully uploaded.
Task i13o8djcoq70bsdut1x73zi69 submitted to the default queue.
Consider tracking the status of the task via CLI: 
        inductiva tasks list --id i13o8djcoq70bsdut1x73zi69
Or, tracking the logs of the task via CLI:
        inductiva logs i13o8djcoq70bsdut1x73zi69
Task i13o8djcoq70bsdut1x73zi69 successfully queued and waiting to be picked-up
for execution...
```
The simulation should take around **2 minutes** to complete and the resulting data will be stored inside a directory within `inductiva-output/{task-id}`. This output data includes some log file, notably `stderr.txt`, `stdout.txt` and `log/SPH_log.txt` detailing the simulation process and any errors that may have occurred. More importantly, the output directory will contain a `vtk` subdirectory where several `vtk` files store data about the fluid particles at each simulation timestep:

```
$ tree inductiva-output/i13o8djcoq70bsdut1x73zi69
inductiva-output/i13o8djcoq70bsdut1x73zi69
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

These `vtk` files are the seeds of our dataset.

## Getting a sense for realistic computational costs
Running the above "base case" simulation was relatively quick since there were only approximately _2,000 particles_ representing the fluid block and the container. However, in the
[study by Sanchez-Gonzalez et al.](https://arxiv.org/abs/2002.09405) that we're
building upon, the researchers modelling 
**significantly larger systems**, involving more than _20,000 particles_. 

> So, what happens when we increase the number of particles in our simulation to a number of that is compatible with the one used by the authors? 

To investigate, let's return to our input directory to revisit the `JSON` file that contains our 
simulation's parameters. To align our simulation more closely with the particle 
count seen in the study, we'll adjust the `particleRadius` from _0.01_ to _0.008_. 
This change will increase the number of particles in our simulation, thus enhancing the overal 
resolution and (hopefully) leading to more physically accurate results.

```{code-block} json
:caption: Set the particle radius to 0.008 in "Configuration" block.
{   
    "Configuration": {
        "stopAt": 4,
        "timeStepSize": 0.01,
        "particleRadius": 0.008,
    }
} 
```

Following the same steps we took to[run our preconfigured "base case"] above, we will see that 
this simulation should take around **30 minutes** to run, and generates roughly _300 MB of data_ 
from _601 frames_. If we're aiming to run 10,000 or even 20,000 variations of our "base case" using the current hardware option, then we will have to wait for many thousands hours to achieve a large-enough dataset. 

To speed this up, we will defintely need to choose much more powerful computational resources to run the required number of simulation. Fortunatelly, this is extremely easy to do using the Inductiva API.


