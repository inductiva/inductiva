---
myst:
  html_meta:
    description: "Learn how to generalize the 'base case' simulation script using Inductiva's Templating Engine to modify parameters programmatically directly from your Python code."
    keywords: "Inductiva API, Programming, HPC, Simulation, Tutorial, Synthetic Data Generation, Physics-ML, SPH"
---

# Generalize the "Base Case" Physical Parameters

In previous chapters of this tutorial, we [set up a "base case" simulation]({% post_url 2024-03-13-api-synthetic-data-generation-2 %}) of a simple "fluid in a box" scenario using SPlisHSPlasH, dropping a 
_0.5m block of water_ into a box from one of its top corners and letting it flow 
for _4 seconds_. Then, we zeroed in on the importance of [choosing the right hardware setup]({% post_url 2024-03-13-api-synthetic-data-generation-3 %}) by testing our "base case" performance 
on different VM types. We compared the beefy but pricey _`c3-standard-88`_ machines 
against the default but more affordable _`c2-standard-4`_ VMs made available by 
our API via Google Cloud, noting the difference in their cost and total execution time.

In this step, we'll expand upon our simple "base case" and **generalize the simulation script**
to be able to programatically set the physical parameters of the simulation.
We will achieve this through using Inductiva's own [Templating Engine](https://docs.inductiva.ai/en/latest/explore_api/templating.html)

## What is Templating?

Templating is a powerful tool that allows you to start with a specific simulation 
file – like our “base case” – containing fixed values for the parameters you wish to 
explore and transform those fixed values into variables that you can now change 
programmatically from your Python code before you submit the simulation for remote 
execution.

The power of templating happens through a simple substitution process. You replace 
the numeric or categorical values in your configuration file with placeholders 
in the following format `{% raw %}{{ variable_name }}{% endraw %}`. This ensures 
that each occurrence of `variable_name` expression within the file is substituted 
with whatever value assigned to it. For example, let's say you embed `{% raw %}variable = {{ value }}{% endraw %}` 
in your template and assign `value = 10`, when you render this template, the placeholder 
is replaced with the assigned value, resulting in `variable = 10` in the final file. 
This process not only generalizes your configuration file but also allows you to 
adjust the values of each variable programmatically via your Python script.

Now, let's revisit our "base case" script to identify the parameters and values 
we want to "generalize". The **first step is to generalize the parameters directly related to the physical properties** of the simulation case itself, like initial conditions, viscosity, or other physical 
description.

## Generalizing our "Base Case" Physical Parameters

As mentioned, we're particularly interested in generalizing the physical properties 
and initial conditions of the fluid block in our simulation. The key parameters 
we want to transform into variables include the fluid block's ***dimensions***, 
***initial position*** and ***initial velocity***. At the same time, we want to 
be able to change the ***density*** and ***viscosity*** of the fluid itself to
to create more diverse examples for our target ML task.

>Let's **<a href="/assets/files/splishsplash-template-dir.zip" download="splishsplash-template-dir.zip" class="bi bi-cloud-download-fill">
<span> download our pre-configured template folder,</span>
      </a>** and store it in a local directory.

If you recall, the configuration for our simulation, including the parameters we're now 
making variable, are located in a [`JSON` file]({% post_url 2024-03-13-api-synthetic-data-generation-2 %}) 
and saved in the local directory as part of the download folder above.

For this dataset, we're keeping the dimensions of the encasing box fixed to focus 
on how changes to the fluid block and its properties impact the simulation 
outcomes without introducing additional variables related to the surrounding 
environment.

Here's an overview of how our templated configuration file looks like, keeping in
mind some sensible defaults we've provided for each variable to ensure the simulation 
can run effectively.

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

After making these changes in our configuration file, executing our simulation 
with different parameters becomes remarkably straightforward. We can now invoke 
the below script and easily fill in the variable values, enabling us to simulate 
a variety of fluid behaviors and transition from modeling water to simulating 
another viscous fluid with a significant horizontal velocity. Here's how:

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

Dynamics of a viscous fluid simulated with **26,600** particles performed using SPlisHSPlasH via the Inductiva API

## Up Next: Adjusting our Hyperparameters

After we've generalized and tested the parameters defining our simulation's 
physical properties, we're now ready to focus on generalizing the hyperparameters 
that affect simulation performance. Key settings, such as simulation resolution 
and iteration count significantly influence both the accuracy of our simulations 
and their associated costs. 

Our [next step]({% post_url 2024-03-24-api-synthetic-data-generation-5 %})
in this tutorial involves an exploration of these hyperparameters, 
not just to generalize and adjust them, but to understand their effect on computational 
cost. This is especially important when considering the potential cost implications 
of running simulations across thousands of configurations. Through our API, we will adjust 
these settings and learn how to reach an acceptable compromise between the fidelity 
of our simulation and the computational cost.
