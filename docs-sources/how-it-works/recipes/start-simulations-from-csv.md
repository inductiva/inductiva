# How to start your simulations from a CSV file

At Inductiva, we aim to allow researchers and engineers to run their simulations
at scale. This can involve, hundreds or even thousands of simulations over a wide
range of parameters. To facilitate this, we provide a simple way to start your
simulations from a CSV file. This allows you to define your simulation parameters
in a structured format, facilitating the process of starting a batch of simulations.

## How to Start Simulations from a CSV File

Here it's a python script that will allow you to start one simulation for each
row in a CSV file:

```python
import inductiva
import csv

# Allocate cloud machine on Google Cloud Platform
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c2d-highcpu-4",
    num_machines=3,
    spot=True)

# Initialize the Simulator
splishsplash = inductiva.simulators.SplishSplash(\
    version="2.13.0")

# Replace 'data.csv' with your actual CSV filename
filename = 'data.csv'

with open(filename, mode='r', newline='') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        density0 = row['density0']
        particleRadius = row['particleRadius']
        viscosity = row['viscosity']

        inductiva.TemplateManager.render_dir( \
            source_dir="my_templates",
            target_dir="my_rendered_files",
            overwrite=True,
            density0=density0,
            particleRadius=particleRadius,
            viscosity=viscosity)

        task = splishsplash.run( \
            input_dir="my_rendered_files",
            sim_config_filename="DamBreakModel.json",
            on=cloud_machine)
        
        # set metadata for the task
        task.set_metadata(
            {
                "density0":f"{density0}",
                "particleRadius":f"{particleRadius}",
                "viscosity":f"{viscosity}"
            }
        )
```

This recipe requires you to use our [templaing system](https://docs.inductiva.ai/en/latest/intro_to_api/templating.html)
to create the template files for your simulation.
