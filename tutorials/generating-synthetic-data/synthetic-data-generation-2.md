---
myst:
  html_meta:
    description: "Dive into the first crucial step of synthetic data generation and learn how to define your 'base case' simulation model."
    keywords: "Inductiva API, Programming, HPC, Simulation, Tutorial, Synthetic Data Generation, Physics-ML, SPH"
---

# Set Up the "Base Case"

In the [introduction]({% post_url 2024-03-13-api-synthetic-data-generation-1 %})
of this tutorial, we outlined and touched on the series of steps we take towards 
generating synthetic datasets for training Physics-ML models. In this chapter, we 
will dive into the first step of this process: **defining a "base case"** simulation 
model of the system we wish to study.

As you recall, we followed through a [study by Sanchez-Gonzalez et al.](https://arxiv.org/abs/2002.09405)
where a Graph Neural Network (GNN) is used to learn and potentially replace classic 
fluid dynamics simulators. The authors of that study were specifically interested in using 
machine learning to mimic Smoothed Particle Hydrodynamics (SPH), a mesh-free fluid 
simulation technique that models both fluids and surface constraints, like walls 
or obstacles, as vast groups of particles that move and interact with each other 
through time and space, according to the laws of fluid dynamics.

The computational challenge of such a simulation method lies in managing 
these interactions, especially with simulations involving millions of particles. 
The study used Graph Neural Networks (GNNs) to model these particles 
and their interactions as a graph, learning to update each particle's acceleration
over time. The researchers use a straightforward simulation scenario of a fluid 
block dropped inside a box, and then model how this fluid interacts with the box's 
surfaces until it stabilizes. Their approach sought to enable GNNs to generalize 
fluid behavior predictions across a range of physical conditions and properties, 
including density, viscosity, and different initial states of position and velocity.

Over this tutorial series, we want to show you how to create a synthetic dataset to train 
this GNN model for such a simulation use case, through the Inductiva API. This all 
begins with preparing our "base case", using one of the open-source Smoothed Particle 
Hydrodynamics (SPH) simulations already built into our API: [SPlisHSPlasH](https://docs.inductiva.ai/en/latest/simulators/SPlisHSPlasH.html). 

<div style="display: flex; justify-content:center">
<video width=500 loop muted autoplay preload="auto">
<source src="../_static/generating-synthetic-data/fluid_drop.mp4" type="video/mp4">
</video>
</div>

Video 1: Dynamics of a water
block simulated with _29791 SPH particles_. Here, the block is left to fall with
gravity - *splash!* Simulation performed via Inductiva API.

## Step 1: Setting Up Our "Base Case" 

Building on the simple simulation scenario used by the researchers in their study, 
our "base case" simulation is quite the same: _a 0.5m cube of water, initially at 
rest at one of the top corners of a sealed 1m cubic box, is dropped at the simulation 
onset, allowing the water to spill and splash against the walls of the closed box 
for 4 seconds._

Make sure you're equipped with the Inductiva API and have your API key 
at hand before diving into this tutorial. If this is your first encounter with our 
API, we highly recommend going through our [Quickstart Tutorial](https://docs.inductiva.ai/en/latest/get_started/installation.html) to set up your environment correctly. 

### Configuring the "Base Case" Parameters

To kick things off, we've pre-configured an input directory containing all 
the configuration files necessary to run the SPlisHSPlasH simulation, and defined 
the particle radius to _0.01_ to decrease the simulation's runtime and test our API. 

>Let's **<a href="/assets/files/splishsplash-base-dir.zip" download="splishsplash-base-dir.zip" class="bi bi-cloud-download-fill">
<span> download our pre-configured input folder,</span>
      </a>** and store it in a local 
      directory.
      
In this folder, we'll find:

- An `.obj` file which is a simple data format that represents the 3D geometry of 
the containers. For our example, it contains a detailed outline of our cubic box
geometry.

- A `JSON` file which contains the simulation parameters required for SPlisHSPlasH to run. 
This script essentially comprises **four key blocks** that define the whole simulation; 
Configuration, RigidBodies, Materials. and FluidModels. This cofiguration file is 
crucial as it provides the flexibility to tweak simulation parameters, enabling 
us to to programmatically run multiple variations of our "base case" and generate 
a diverse enough dataset to train a machine learning (ML) model.

Let's take a closer look at these parameters and what they mean:

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
It's essential to understand that these parameters collectively 
control the simulation's behavior. Most importantly, the **particle size plays 
a crucial role in determining the computational requirements of the simulation.**

### Running our "Base Case"

To run our "base case", we'll initialize the SPlisHSPlasH simulator using the API 
with just a couple lines of code:

```python
import inductiva

# Set path to the input directory with the SPlisHSPlasH configuration files
input_dir = "splishsplash-base-dir"

# Initialize the SPlisHSPlasH simulator through the API
splishsplash = inductiva.simulators.SplishSplash()

# Run the simulation task with the parameters defined in the .json file
task = splishsplash.run(input_dir=input_dir,
                        sim_config_filename="config.json")

# Wait for the simulation to complete and download the outputs
task.wait()
task.download_outputs()
```
When launching the simulation, the API will upload the input data from our local 
directory and schedule the simulation `task` for execution. It will also provide 
details about the `task`, including its ID and the machine group assigned for 
computation:

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
Simulation metadata logged to: inductiva_output/task_metadata.json
Task i13o8djcoq70bsdut1x73zi69 configurations metadata saved to the tasks
metadata file task_metadata.json in the current working directory.
Consider tracking the status of the task via CLI: 
        inductiva tasks list --task-id i13o8djcoq70bsdut1x73zi69
Or, tracking the logs of the task via CLI:
        inductiva logs i13o8djcoq70bsdut1x73zi69
Task i13o8djcoq70bsdut1x73zi69 successfully queued and waiting to be picked-up
for execution...
```
The simulation should take around **2 minutes** to complete and the resulting 
data will be stored inside a directory within `inductiva-output/{task-id}`. This 
includes `stderr.txt`, `stdout.txt` and `log/SPH_log.txt` log files detailing the 
simulation process and any encountered errors, along with `vtk` files stored in 
the `vtk` directory that contain particle data at each simulation timestep:

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

Running the above "base case" simulation was relatively quick due to its low 
resolution count of approximately _2,000 particles_. However, in the
[study by Sanchez-Gonzalez et al.](https://arxiv.org/abs/2002.09405) that we're
building upon, the researchers employ a Graph Neural Network (GNN) to simulate 
**significantly larger systems**, involving more than _20,000 particles_. 

### Enhancing the "Base Case" Resolution

Let's return to our input directory to revisit the `JSON` file that contains our 
simulation's parameters. To align our simulation more closely with the particle 
count seen in the study, we'll adjust the `particleRadius` from _0.01_ to _0.008_. 
This change will increase the number of particles in our simulation, enhancing the 
simulation's resolution for more detailed results.

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

Following the same steps we took to [run our preconfigured "base case"](#running-our-base-case) above, 
this simulation should take around **30 minutes** to run, generating roughly _300 MB of data_ 
from _601 frames_. While tweaking our "base case" to increase particle count gave 
us more detailed data, it also meant waiting 30 minutes for a single run simulating 
4 seconds of fluid dynamics. If we're aiming to run 10,000 or even 20,000
variations of our "base case" to have a large enough dataset, we're looking 
at an astronomical number of computation hours! To speed this up, choosing more 
powerful computational resources becomes an essential step.

## Up Next: Choosing more Powerful Hardware

In this step, we've learned how to set up our "base case" simulation and configure 
its parameters. Next, we'll dive into how we can speed things up and keep the waiting 
to a minimum, particularly as we tweak the parameters to simulate a much higher 
particle count. Through our API, we'll explore available resources, their cost, 
and learn how to choose the right machine setup for our needs.

See you in the [next one]({% post_url 2024-03-13-api-synthetic-data-generation-3 %})!

