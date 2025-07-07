# Creating a SPlisHSPlasH Dataset

In the previous sections, we introduced our templating system and showed how
easy it is to run multiple simulations in parallel to explore a wide range of
parameters.

In this part of the tutorial, we’ll put that into practice.

<p align="center"><img src="../../_static/combined.gif" alt="Visualization of 10 simulations" width="600"></p>

We will run 400 simulations while varying the following parameters:

* **Particle start and end positions**: randomly shifted within ±0.33 of the original position, moving the particles up/down
* **Initial velocity**: randomly set between -1 and 1 in all axes
* **Density**: randomly selected between 800 and 1200
* **Viscosity**: randomly chosen between 0.05 and 1.05

Each parameter will be sampled randomly within the specified bounds, allowing us to generate a diverse and representative dataset.

## The code

In order to run this dataset generation, we will use the following code:

```python
import inductiva
import random

# Allocate cloud machine on Google Cloud Platform
cloud_machine = inductiva.resources.ElasticMachineGroup(
    provider="GCP",
    machine_type="c2d-highcpu-16",
    min_machines=1,
    max_machines=50)

# Initialize the Simulator
splishsplash = inductiva.simulators.SplishSplash()

# Assuming the template folder was downloaded to the local directory,
# set the path to it
template_dir = "splishsplash-base-dir"

# Parameters to explore
start=[]
end=[]
initialVelocity=[]
density0=[]
viscosity=[]

for i in range(400):
    # random scale where the volume stays as 0.5*0.5*0.5
    
    # random number between -0.33 and 0.33 for y-axis
    rand_y = ((random.random() * 2 ) - 1) / 3

    start.append([-0.5, -0.5+rand_y, -0.5])
    end.append([0.5, 0.5+rand_y, 0.5])
    
    # random initial velocity between -1 and 1
    rand_s_x = ((random.random() * 2 ) - 1) / 1
    rand_s_y = ((random.random() * 2 ) - 1) / 1
    rand_s_z = ((random.random() * 2 ) - 1) / 1
    initialVelocity.append([rand_s_x,rand_s_y,rand_s_z])

    # random density between 800 and 1200
    rand_d= ((random.random()*2) -1 )*200
    density0.append(1000 + rand_d)

    
    rand_viscosity = (random.random() * 1)
    # random viscosity between 0.05 and 1.05
    viscosity.append(0.05+rand_viscosity)
    

for s,e,v,ini_vel,d in zip(start,end,viscosity,initialVelocity,density0):
    # Define the directory where the rendered templates will appear filled 
    # with the values of the variables defined below
    target_dir = "splishsplash-hyperparameter-search"
    inductiva.TemplateManager.render_dir(
        source_dir=template_dir,
        target_dir=target_dir,
        particle_radius=0.015,
        start=s,
        end=e,
        initial_velocity=ini_vel,
        density=d,
        viscosity=v,
        # overwrite=True to allow overwriting the same directory
        overwrite=True)

    task = splishsplash.run(
        input_dir=target_dir,
        sim_config_filename="config.json",
        project="splishsplash_400",
        on=cloud_machine,
        resubmit_on_preemption=True)
    
    # Set metadata for the task
    task.set_metadata(
        {
            "start": str(s),
            "end": str(e),
            "initial_velocity": str(ini_vel),
            "density0": str(d),
            "viscosity": str(v),
            "particle_radius": str(0.015),

        }
    )
```

We start this script by allocating a cloud Elastic Machine Group. This resource
will have a minimum of 1 and a maximum of 50 machines, this will allow us to run
up to 50 simulations in parallel. Once the simulations start to finish our Machine
Group will scale down, keeping the costs to a minimum.

Following that we have a loop that will generate the random parameter values for
the whole 400 simulations.

Lastly we have another loop that will start each simulation, using our templating
system to replace the parameters in the simulation files with the generated values
and submitting the simulation to the cloud machine. We also added some metadata
to each simulation just to keep track of the parameters used in each simulation.

This is all you need to generate a dataset with 400 simulations, each with
different parameters.

## Results

Once all the simulation are finished, you can download the results from all the simulations
by running the following terminal command:

```bash
inductiva project download splishsplash_400
```

You can also run this script to get some information about your project:

```
import inductiva

project = inductiva.projects.Project("splishsplash_400")

print(project)

Project 'splishsplash_400' created at 2025-07-04 13:31.

Total number of tasks: 400

Tasks by status:
  TaskStatusCode.SUCCESS: 400

Estimated total computation cost: 2.55 US$
```

From this output we can see that our project has 400 tasks, all of them
successfully completed, and the estimated total computation cost is around 2.55 US$.

Keep following the tutorial to learn how you can visualize the results of these simulations!