# Run Multiple Simulations in Parallel

Running multiple simulations in parallel can significantly reduce waiting times, 
especially useful when exploring various parameter values or running a large number 
of simulations for a sensitivity analysis. This how-to guide will walk you through 
using Machine Groups to run several simulations in parallel, using the 
[templating mechanism](https://tutorials.inductiva.ai/intro_to_api/templating.html)
integrated within the Inductiva API. 
This approach makes it easy to explore variations of a base simulation scenario. 
As a practical example, we will use a coastal dynamics simulation with 
the [SWASH simulator](https://tutorials.inductiva.ai/simulators/SWASH.html).

## 1. Setting Up Your Machine Group

First, create a MachineGroup to run your simulations in parallel:

```python
import inductiva

# Instantiate a MachineGroup object with 5 preemptible machines of type
# c2-standard-30 and start it immediately
machine_group = inductiva.resources.MachineGroup(
    machine_type="c2-standard-30", num_machines=5, spot=True)
machine_group.start()
```

## 2. Preparing Simulation Inputs

Download and prepare the input files for your simulations:

```python
# Download input files for the SWASH simulation
template_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/swash-template-example.zip", unzip=True)
```
## 3. Running the Simulations

Define the variations for your simulation - _here, different water levels_ — and 
launch the simulations:

```python
# Initialize the SWASH simulator
swash = inductiva.simulators.SWASH()

# Define different water levels to explore
water_levels_list = [3.5, 3.75, 4.0, 4.5, 5.0]

# Launch multiple simulations for each water level
for i, water_level in enumerate(water_levels_list):
    target_dir = f"./inductiva_input/swash-sim-{i}"  
    inductiva.TemplateManager.render_dir(
                            source_dir=template_dir,
                            target_dir=target_dir,
                            water_level=water_level,
                            overwrite=False)

    # Run the simulation on the dedicated MachineGroup
    task = swash.run(input_dir=target_dir,
                    sim_config_filename="input.sws",
                    on=machine_group)
```

## 4. Monitoring Simulations

The template mechanism will allow you to explore 5 different variations of the
simulation, each with a different water level. The simulations will be submitted
to our dedicated machine group and will run in parallel.

You can check the status of these simulations through the
[Inductiva CLI](../cli/streaming-logs.md), and you'll see that it took only **1 minute** 
from the moment they were submitted until they start running:

```bash
$ inductiva tasks list

    ID                         Simulator    Status    Submitted         Started           Computation Time    Resource Type

    57mr4kas99jxb9titkeackano  swash        started   01 Feb, 09:07:19  01 Feb, 09:08:03  *0:03:12            c2-standard-30
    ox8718m0pwfi02zczui3qky4w  swash        started   01 Feb, 09:07:17  01 Feb, 09:08:02  *0:03:14            c2-standard-30
    mak1ji62s7axf7mespkc36g7e  swash        started   01 Feb, 09:07:15  01 Feb, 09:08:03  *0:03:14            c2-standard-30
    ijyu8bkvme7vg9k0kj6v23gxa  swash        started   01 Feb, 09:07:14  01 Feb, 09:08:02  *0:03:16            c2-standard-30
    g5qq5c9mk2nr5wqhzef38sdm4  swash        started   01 Feb, 09:07:12  01 Feb, 009:08:01  *0:03:17            c2-standard-30
```

Running simulations in parallel significantly optimizes the use of computational 
resources. In our example, all five simulations start in parallel and complete in 
**9 minutes and 55 seconds**, roughly the time it would take to run one single 
simulation!

````{eval-rst}
.. important::
   After completing the simulations, remember to release your computational resources 
   to avoid unnecessary charges!
````
```bash
# Terminate all computational resources
$ inductiva resources terminate --all
```

