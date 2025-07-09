# Generalizing the Base Case
In this step, we’ll build on our simple base case and generalize the simulation script to allow programmatic control over the physical parameters of the simulation. To do this, we’ll use Inductiva’s built-in [Templating System](https://website-staging.inductiva.ai/guides/documentation/intro/templating).

## What is Templating?
Templating is a powerful mechanism that lets you start with a specific simulation file, such as our base case, that contains fixed values for the parameters you want to explore. You can then transform those fixed values into variables that can be modified directly from your Python code before submitting the simulation for remote execution.

Now, let’s revisit our base case script and identify the parameters and values we want to generalize. We'll start by focusing on those directly related to the physical properties of the simulation, such as initial conditions, viscosity, and other aspects of the fluid’s physical behavior.

## Generalizing Physical Parameters
The key parameters we want to convert into variables include the fluid **block’s dimensions**, **initial position**, and **initial velocity**. Additionally, we want to make the **fluid's density** and **viscosity** configurable in order to generate more diverse examples for our target machine learning task.

Start by [downloading our pre-configured template folder](https://storage.googleapis.com/inductiva-api-demo-files/splishsplash-template-dir.zip) and saving it to a local directory. As you may recall, the simulation configuration, including the parameters we are now making variable, is stored in a `.json` file, which is included in the downloaded folder.

Below is an overview of our templated configuration file. Note that we’ve defined sensible default values for each variable to ensure the simulation runs smoothly, even before custom values are provided.

```text
{
    "Configuration": {
        "stopAt": 4,
        "timeStepSize": 0.01,
        "particleRadius": 0.008,
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
    "RigidBodies": [
        {
            "geometryFile": "unit_box.obj",
            "translation": [0, 0, 0],
            "scale": [1, 1, 1],
            "isDynamic": false
        }
    ],
    "Materials": [
        {
            "id": "Fluid",
            "density0": {{ density | default(1000) }},
            "viscosity": {{ viscosity | default(1e-6) }},
            "viscosityMethod": 6
        }
    ],
    "FluidModels": [
        {
            "id": "Fluid",
            "particleFile": "unit_box.obj",
            "translation": {{ initial_position | default([0, 0, 0]) }},
            "scale": {{ dimensions | default([0.5, 0.5, 0.5]) }},
            "initialVelocity": {{ initial_velocity | default([0, 0, 0]) }}
        }
    ]
}
```

Executing our simulation with different parameters is remarkably straightforward. All it takes is using Inductiva’s `TemplateManager` class. In the script below, you can see how easily we populate the template variables by calling the `TemplateManager.render_dir()` method:

```python
import inductiva

# Allocate cloud machine on Google Cloud Platform
cloud_machine = inductiva.resources.MachineGroup(
    provider="GCP",
    machine_type="c2d-highcpu-4",
    spot=True)

# Set the path to the local directory where the template folder was downloaded
template_dir = "splishsplash-template-dir"

# Specify the initial conditions for the fluid simulation and properties
# Example: A high horizontal initial speed m/s in the x direction 
initial_velocity = [4, 0, 0]
# Represents a fluid with higher viscosity m^2/s
kinematic_viscosity = 2       
# A denser fluid compared to water 1000 kg/m^3
density = 2500                

# Define the directory where the rendered templates will appear filled 
# with the values of the variables defined above
rendered_dir = "./splishsplash-viscous-fluid/"

inductiva.TemplateManager.render_dir(
    source_dir=template_dir,
    target_dir=rendered_dir,
    density=density,
    viscosity=kinematic_viscosity,
    initial_velocity=initial_velocity,
    overwrite=True)

# Initialize the Simulator and run the simulation
SPlisHSPlasH = inductiva.simulators.SplishSplash()
task = SPlisHSPlasH.run(
    input_dir=rendered_dir,
    sim_config_filename="config.json",
    on=cloud_machine)

# Wait for the simulation task to complete and download the results
task.wait()

# Ensure that the allocated resources are terminated
# This is crucial to avoid incurring unnecessary costs from lingering resources
cloud_machine.terminate()

task.download_outputs()
```

And that’s it! Inductiva’s templating mechanism allows us to simulate a wide range of fluid behaviors, from modeling water to 
simulating viscous fluids with horizontal velocity, using essentially the same configuration file. It really is as simple as it seems.

Now, we’re ready to take things a step further by generalizing the simulator’s hyperparameters, which have a direct impact on the 
execution time of each simulation. Our focus will be on **particle radius**, a key setting that determines the total number of particles used.

This hyperparameter plays a crucial role, as it significantly affects both the accuracy of the results and the overall runtime.