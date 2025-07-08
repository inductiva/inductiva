# Generate the Dataset
In the previous sections, we built a base SPH simulation using **SPlisHSPlasH**, used the template system to configure simulation parameters to reflect real-world variability, and demonstrated how easy it is to run multiple simulations in parallel to explore a wide range of parameters.

In this final section, we’ll put all of that into practice and show you how to scale up to generate synthetic data in bulk.

<p align="center"><img src="../../_static/combined.gif" alt="Visualization of 10 simulations" width="600"></p>

We will run **400 simulations**, randomly varying the following parameters:
* **Particle start and end positions**: randomly shifted within ±0.33 of the original position, moving the particles up/down
* **Initial velocity**: randomly set between -1 and 1 in all axes
* **Density**: randomly selected between 800 and 1200
* **Viscosity**: randomly chosen between 0.05 and 1.05

Each parameter will be sampled independently within its specified range, allowing us to create a diverse and representative synthetic dataset.

## Code Overview
Generating the dataset is as simple as executing the Python script shown below:

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

We begin the script by allocating a cloud [Elastic Machine Group](https://inductiva.ai/guides/scale-up/parallel-simulations/set-up-elastic-machine-group) with a minimum of 1 and a maximum of 50 machines. This setup allows 
us to run up to 50 simulations in parallel. As simulations complete, the Machine Group automatically scales down, helping to keep costs to 
a minimum.

Next, a loop generates random parameter values for all 400 simulations.

Lastly, another loop starts each simulation by using our templating system to replace the parameters in the simulation files with the generated values and submits the simulation to the cloud machine. Metadata is also added to each simulation to keep track of the parameters used.

And that’s it — with this setup, you can generate a dataset of 400 simulations, each with its own unique parameter configuration.

## Results
Once all the simulations have finished, you can download the results from all the simulations by executing the following terminal command:

```bash
inductiva project download splishsplash_400
```

You can also run this script to get detailed information about your project:

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

From the output, we can see that the project contains **400 tasks**, all **successfully completed**, with an estimated total computation cost of **2.55 US$**.

(remeter para os tutoriais de viz)
