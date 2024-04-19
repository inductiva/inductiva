---
myst:
  html_meta:
    description: "Explore the impact of adjusting simulation hyperparameters like particle radius on the simulation's computational cost and runtimes."
    keywords: "Inductiva API, Programming, HPC, Simulation, Tutorial, Synthetic Data Generation, Physics-ML, SPH"
---

# Test the Impact of Hyperparameter Changes

In the preceding step of our tutorial, we used Inductiva's [Templating Engine](https://docs.inductiva.ai/en/latest/explore_api/templating.html).
to generalize the physical properties and initial conditions of the fluid block
in our "base case" simulation. By substituting certain values in our configuration
file with placeholders, we enabled the template to incorporate assigned values,
enabling programmable adjustments of each variable via our Python script.

Now, we will move on to the hyperparameters of our simulation, which are important
due to their influence on both the fidelity of the simulation and the associated
computational requirements and costs. Among various hyperparameters, we will specifically
focus on adjusting the ***particle radius*** value and how it affects computational
costs and data output. To illustrate this, we will run four simulations with decreasing
particle radii while keeping the remaining parameters fixed to understand
these impacts better.

## Generalizing the `Particle Radius`

Using the same methodology from our previous step, we will substitute the categorical
and numerical value of the `particle radius` hyperparameter in the `.JSON` configuration
file we previously downloaded with a templated variable:

_`{% raw %}{{ variable_name }}{% endraw %}`_

Here's what our templated configuration file looks like, where we've kept all
hyperparameters fixed except for the particle radius:

```text
"Configuration": {
    "stopAt": 1,
    "timeStepSize": 0.001,
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

Now, let's save the `.JSON` file in the local directory, within the download folder,
to prepare for running the simulation.

## Running the Templated Simulation

After updating our configuration file, we're ready to run simulations with different
particle sizes. We will invoke the below Python script to set up a group of machines
for the simulations and initialize our templating engine to easily fill in the
`particle size` variable values. We're testing four sizes: **0.01, 0.008, 0.006, and 0.004 meters.**

```python

import inductiva

# Launch a machine group with a c3-standard-4
machine_group = inductiva.resources.MachineGroup("c3-standard-4")
machine_group.start()

# Assuming the template folder was downloaded to the local directory,
# set the path to it.
template_dir = "./splishsplash-template-dir"

# Define the radii of the particles
particle_radii = [0.01, 0.008, 0.006, 0.004]

tasks_list = []

# Initialize the templating engine
template_manager = inductiva.TemplateManager(template_dir)
for n, radius in enumerate(particle_radii, start=1):
    # set the output directory to a different folder for each simulation
    # and render the entire content of the template directory with the
    # particle radius set to the current value
    template_manager.set_root_dir("splishsplash-hyperparameter-search_%d" % n)
    template_manager.render_dir(particle_radius=radius)

    task = SPlisHSPlasH.run(input_dir=template_manager.get_root_dir(),
                            sim_config_filename="config.json",
                            on=machine_group)
    tasks_list.append(task)

```

For each particle size, our API creates a new folder, updates the settings with
the new size, and starts the simulation. This lets us run all four simulations at
the same time, making it faster to see how changing the particle size affects the
results.

Instead of waiting for each simulation to complete within the running script,
we can track their progress and collect the results afterwards using our [Command Line Interface (CLI)](https://docs.inductiva.ai/en/latest/cli/cli-overview.html):

```bash
# Monitor the status of the tasks every 10 seconds
$ inductiva task list -n 4 --watch 10

# Download the results of the tasks
$ inductiva storage list -m 4

```

## Impact of Particle Radius on Simulation Volume and Data Output

The table below shows how the different particle radius values we've chosen affect
the total number of particles and the amount of data generated. As we decrease
the particle radius, we need more particles for the simulation, increasing the
data volume. Specifically, halving the particle radius results in an eightfold
increase in the number of particles and the data produced, due to the volume
scaling with the cube of the particle radius ($$r_{particle}{^3}$$).

<span class="mt-0 block sm:text-left text-base"><strong>Table 1.</strong> Number
of particles and total size of simulation output for each particle radius.</span>

| Particle Radius | Total N. of Particles | Size of data produced  |
|-----------------|-----------------------|------------------------|
| 0.01            | 15625                 | 166 MB                 |
| 0.008           | 29791                 | 213 MB                 |
| 0.006           | 73720                 | 525 MB                 |
| 0.004           | 244584                | 1.76 GB                |

## Performance and Cost Analysis Across Machine Configurations

Reducing the particle radius naturally results in a higher number of particles
required for the simulation to run, and this demands more computing power.
While upgrading to more powerful machines can reduce computational times, it also
incurs significantly higher costs. To better understand these factors, we've run
our simulations across various machine configurations. The resulting data, shown
in the tables below, outlines runtime and costs for our four simulations
on three different machine types:

- `c3-standard-4` at **$0.23 per hour**
- `c3-standard-8` at **$0.459 per hour**
- `c3-standard-88` at **$5.053 per hour**

<span class="mt-0 block sm:text-left text-base"><strong>Table 2.</strong>
Simulation runtimes and costs on a `c3-standard-4` machine, priced at $0.23 per
hour, as performed via the Inductiva API.</span>

| Particle Radius | Time to run | Cost in $      |
|-----------------|-------------|----------------|
| 0.01            |   16m31s    | 0.06           |
| 0.008           |   27m11s    | 0.10           |
| 0.006           |   56m59s    | 0.22           |
| 0.004           | 6h11m27s    | 1.42           |

<span class="mt-0 block sm:text-left text-base"><strong>Table 3.</strong>
Simulation runtimes and costs on a `c3-standard-8` machine, priced at $0.459 per
hour, as performed via the Inductiva API.</span>

| Particle Radius | Time to run | Cost in $      |
|-----------------|-------------|----------------|
| 0.01            |    9m38s    | 0.07           |
| 0.008           |   15m11s    | 0.12           |
| 0.006           |   35m34s    | 0.27           |
| 0.004           | 3h28m31s    | 1.60           |

<span class="mt-0 block sm:text-left text-base"><strong>Table 4.</strong>
Simulation runtimes and costs on a `c3-standard-88` machine, priced at $5.053 per
hour, as performed via the Inductiva API.</span>

| Particle Radius | Time to run | Cost in $      |
|-----------------|-------------|----------------|
| 0.01            |    3m27s    | 0.29           |
| 0.008           |    5m16s    | 0.44           |
| 0.006           |    8m20s    | 0.70           |
| 0.004           | 1h25m04s    | 7.16           |

Notice that choosing more powerful hardware does not correspondingly decrease
computational times as much as it increases costs. Upgrading the virtual CPUs
from 8 to 88 doesn't proportionally reduce computational time 11-fold. **On average, run times typically drop by only about a factor of 3, while costs soar by 5 times.**

## Up Next: Benchmarking Computational Resources

In this step, we assessed how adjusting hyperparameters, particularly the `particle radius`,
impacts the computational times and costs of our "base case" simulation. We ran
four simulations with increasingly smaller particle radii while keeping all other
parameters fixed. This change naturally increased the number of particles needed
to run the simulation and the amount of data it produced, which in turn demanded
more computing power. Following this, we looked into how these changes influenced
the cost and duration of computing by running the simulation across different
hardware configurations.

Performing such benchmarks becomes crucial for understanding the necessary trade-offs
in synthetic data generation, especially considering the potential to scale up to
10,000 simulations in a manner that is both time-efficient and cost-effective.
If we revisit the [study](https://arxiv.org/abs/2002.09405) we built upon for
this tutorial, a higher particle count results in a larger graph, potentially
complicating memory requirements for training graph neural networks (GNNs).
In our next chapter, we'll dive deeper into the computational costs and runtimes
for such large scale-ups and learn how to efficiently benchmark computational resources.
