# Testing the Impact of Hyperparameters
Next, we'll take a closer look at one of the most important hyperparameters in our simulation: **particle radius**. This parameter plays a key role in balancing the fidelity of the SPH simulation with its computational cost. Our goal is to identify a value that produces a particle count comparable to that used by [Sánchez-González et al.](https://arxiv.org/abs/2002.09405), while keeping the computational requirements manageable.

To approach this systematically, we’ll once again use Inductiva’s templating mechanism. This allows us to replace the fixed numerical value of the particle radius in the .json configuration file with a templated variable, making it easy to experiment with different values programmatically.

To perform a hyperparameter search over the particle radius while keeping all other simulation parameters fixed, modify the templated configuration file as shown below:

```text
"Configuration": {
    "stopAt": 4,
    "timeStepSize": 0.01,
    "particleRadius": {{ particle_radius | default(0.008) }},
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
}
```

This setup uses templating to make `particleRadius` the only configurable parameter, with a default value of 0.008.

Now, save the .json file in the local directory, specifically within the download folder, in preparation for running three **parallel** simulations with particle radii of 0.015, 0.012, and 0.010 meters. The following code demonstrates how to set up and execute these simulations:

```python
import inductiva

# Allocate cloud machine on Google Cloud Platform
cloud_machine = inductiva.resources.MachineGroup(
    provider="GCP",
    machine_type="c2d-highcpu-4",
    num_machines=3)

# Initialize the Simulator
splishsplash = inductiva.simulators.SplishSplash()

# Assuming the template folder was downloaded to the local directory,
# set the path to it
template_dir = "/Path/to/splishsplash-template-dir"

# Define the radii of the particles
particle_radii = [0.015, 0.012, 0.010]

for n, radius in enumerate(particle_radii, start=1):
    # Define the directory where the rendered templates will appear filled 
    # with the values of the variables defined below
    target_dir = f"splishsplash-hyperparameter-search_{n}"
    inductiva.TemplateManager.render_dir(
        source_dir=template_dir,
        target_dir=target_dir,
        particle_radius=radius,
        overwrite=False)

    task = splishsplash.run(
        input_dir=target_dir,
        sim_config_filename="config.json",
        project=f"SplishSplash-radius-study",
        on=cloud_machine)
```

For each particle size, our API creates a new folder, updates the configuration with the corresponding radius, and launches the 
simulation. This allows us to run all three simulations in parallel on separate `c2d-highcpu-4` machines, speeding up the process 
of evaluating how particle size impacts the results.

While the simulations are running, you can monitor their progress and download the results afterward using our Command Line Interface (CLI):

```
# Monitor the status of the tasks every 10 seconds
$ inductiva tasks list -n 4 --watch 10

# Download the results of the tasks
$ inductiva projects download SplishSplash-radius-study
```

## Results 

### Particle Radius and Data Generation
The table below illustrates how varying the particle radius affects both the total number of particles and the amount of data generated during a 6-second simulation. As the particle radius decreases, more particles are required to occupy the same volume, leading to a corresponding increase in data output. Notably, halving the particle radius results in an eightfold increase in the number of particles, as expected based on volume scaling.

| Particle Radius | Total nº of Particles | Data Produced |
| --------------- | --------------------- | ------------- |
| 0.015           | 32768                 | 204 MB        |
| 0.012           | 68921                 | 422 MB        |
| 0.009           | 166375                | 1.01 GB       |

For reference, the dataset produced by [Sánchez-González et al.](https://arxiv.org/abs/2002.09405) was based on simulations containing approximately 8,000 to 25,000 particles. This suggests that, to create a comparable dataset, it is likely unnecessary to use a particle radius smaller than 0.015. Additionally, since each simulation can generate hundreds of megabytes of data, reducing the particle radius would substantially increase storage requirements — making it difficult to manage data from thousands of simulations. This also adds complexity to training GNNs on such a large volume of data.

### Runtime and Cost per Particle Radius
Reducing the particle radius naturally increases the number of particles required for the simulation and, consequently, the overall runtime.

The table below shows the runtime and corresponding cost of running simulations at the three particle radii under consideration.

| Particle Radius (m) | Machine Type    | Execution Time | Estimated Cost (USD) |
|---------------------|-----------------|----------------|----------------------|
| 0.015               | c2d-highcpu-16  | 2 min, 57s     | 0.0045               | 
| 0.012               | c2d-highcpu-16  | 6 min, 55s     | 0.011                | 
| 0.009               | c2d-highcpu-16  | 30 min, 40s    | 0.047                | 

As expected, reducing the particle radius, and consequently increasing the number of particles, leads to a significant rise in 
computation time.

Based on performance and cost considerations, committing to a particle radius of 0.015 is reasonable for running over **10,000** 
simulations **in parallel**. At this scale, the total cost would be approximately **US$45**.

This positions **Inductiva** as a fast and cost-effective platform for large-scale Physics-Informed Machine Learning (PIML) model training.