---
myst:
  html_meta:
    description: "Learn how to generalize the 'base case' simulation script using Inductiva's Templating Engine to modify parameters programmatically directly from your Python code."
    keywords: "Inductiva API, Programming, HPC, Simulation, Tutorial, Synthetic Data Generation, Physics-ML, SPH"
---

# Generalizing the "Base Case"

In this step, we'll expand upon our simple "base case" and **generalize the simulation script**
to be able to programatically set the physical parameters of the simulation.
We will achieve this through using Inductiva's own [Templating Engine](https://docs.inductiva.ai/en/latest/explore_api/templating.html)

## What is Templating?

Templating is a powerful tool that allows you to start with a specific simulation  file – like our “base case” – containing fixed values for the parameters you wish to explore and transform those fixed values into variables that you can now change 
programmatically from your Python code before you submit the simulation for remote execution.

The power of templating happens through a simple substitution process. You replace the numeric or categorical values in your configuration file with placeholders in the following format `{% raw %}{{ variable_name }}{% endraw %}`. This ensures that each occurrence of `variable_name` expression within the file is substituted with whatever value assigned to it. 

For example, let's say you embed `{% raw %}variable = {{ value }}{% endraw %}` in your template and assign `value = 10`, when you render this template, the placeholder is replaced with the assigned value, resulting in `variable = 10` in the final configutation file. This process not only generalizes your configuration file but also allows you to adjust the values of each variable programmatically via your Python script.

Now, let's revisit our "base case" script to identify the parameters and values 
we want to "generalize". The **first step is to generalize the parameters directly related to the physical properties** of the simulation case itself, like initial conditions, viscosity, or other physical 
description.

## Generalizing Physical Parameters

The key parameters we want to transform into variables include the fluid block's ***dimensions***, ***initial position*** and ***initial velocity***. At the same time, we want to be able to change the ***density*** and ***viscosity*** of the fluid itself toto create more diverse examples for our target ML task.

>Let's **<a href="/assets/files/splishsplash-template-dir.zip" download="splishsplash-template-dir.zip" class="bi bi-cloud-download-fill">download our pre-configured template folder,</a>** and store it in a local directory.

If you recall, the configuration for our simulation, including the parameters we're now 
making variable, are located in a [`JSON` file](synthetic-data-generation-2.md) and saved in the local directory as part of the download folder above.

For this dataset, we're keeping the dimensions of the encasing box fixed to focus on how changes to the fluid block and its properties impact the simulation outcomes without introducing additional variables related to the surrounding 
environment.

Here's an overview of how our templated configuration file looks like, keeping in mind some sensible defaults we've provided for each variable to ensure the simulation can run effectively.

```text
{
    "Configuration": {
        "stopAt": 4,
        "timeStepSize": 0.001,
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
            "translation": {{ initial_position | default([0.05, 0.05, 0.45]) }},
            "scale": { dimensions | default([0.5, 0.5, 0.5]) },
            "initialVelocity": { initial_velocity | default([0, 0, 0]) }
        }
    ]
}
```

After making these changes in our configuration file, executing our simulation with different parameters becomes remarkably straightforward. We can now invoke the below script and easily fill in the variable values, enabling us to simulate a variety of fluid behaviors and transition from modeling water to simulating another viscous fluid with a significant horizontal velocity. Here's how:

```python
import inductiva

# Launch a machine group with a c3-standard-4
machine_group = inductiva.resources.MachineGroup("c3-standard-4")
machine_group.start()

# Set the path to the local directory where the template folder was downloaded
template_dir = "./splishsplash-template-dir"

# Specify the initial conditions for the fluid simulation and the fluid's properties
initial_velocity = [4, 0, 0]  # Example: A high horizontal initial speed m/s in each direction 
kinematic_viscosity = 2       # Represents a fluid with higher viscosity m^2/s
density = 2500                # A denser fluid compared to water kg/m^3


# Initialize the templating manager, define the root directory name for the
# rendered files and render all template files inside the template directory
# with the values of the variables defined above.
template_manager = inductiva.TemplateManager(template_dir)
template_manager.set_root_dir("splishsplash-viscous-fluid")
template_manager.render_dir(density=density,
                            viscosity=kinematic_viscosity,
                            initial_velocity=initial_velocity)

# Initialize the simulator and run the simulation
SPlisHSPlasH = inductiva.simulators.SplishSplash()
task = SPlisHSPlasH.run(input_dir=template_manager.get_root_dir(),
                        sim_config_filename="config.json",
                        on=machine_group)

# Wait for the simulation task to complete and download the results
task.wait()
task.download_outputs()

# Ensure that the allocated resources are terminated
    # This is crucial to avoid incurring unnecessary costs from lingering resources
machine_group.terminate()
```

The result looks something like this:

<div style="display: flex; justify-content:center">
<video width=500 loop muted autoplay preload="auto">
<source src="../_static/generating-synthetic-data/viscous_flow.mp4" type="video/mp4">
</video>
</div>

Video 2. Dynamics of a viscous fluid simulated with **26,600** particles performed using SPlisHSPlasH via the Inductiva API. 

This is it! It's really as simple as it seems! This is the power of the templating mechanism provided by Inductiva API.

And now we are ready to go one step further to generalize the *hyperparameters* of the simulator, which directly affect the computational cost of the simulation. We will be focusing on two key hyperparamters: i) the size of the particles (which controls the number of particles used in the simulation) and ii) the iteration time step. These two hyperparameters impact significantly both the accuracy of our simulations and their associated costs.