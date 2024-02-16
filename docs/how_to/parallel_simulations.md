# Run Simulations in Parallel

Machine groups can be used to run multiple simulations in 
parallel, by having each of the machines in the group run a separate simulation.
This is a powerful mechanism to reduce waiting times when you need to run a large
number of simulations, such as when you are exploring different parameter values
or running a large number of simulations for a sensitivity analysis. 

To exemplify, we will use the [templating mechanism](../introduction/templating.md)
that is built in the Inductiva API to explore variations of a base simulation scenario.
More specifically, we will consider a coastal dynamic simulation using the SWASH
simulator, where we want to explore the effect of different water levels on the
simulation results. 

As an example, let's run 5 different simulations in parallel. 

```python
import inductiva
from inductiva import mixins

# Instantiate a MachineGroup object with 5 preemptible machine of type
# c2-standard-30 and start it immediately
machine_group = inductiva.resources.MachineGroup(
    machine_type="c2-standard-30", num_machines=5, spot=True)
machine_group.start()

# Download the input files for the SWASH simulation
input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "swash-template-example.zip", unzip=True)

# Initialize the template file manager
file_manager = mixins.FileManager()

# Initialize the SWASH simulator
swash = inductiva.simulators.SWASH()

# Explore the simulation for different water levels
water_levels_list = [3.5, 3.75, 4.0, 4.5, 5.0]

# Launch multiple simulations
for water_level in water_levels_list:
    # Set the root directory and render the template files into it.
    file_manager.set_root_dir("swash-input-example")
    file_manager.add_dir(input_dir, water_level=water_level)

    # Run the simulation on the dedicated MachineGroup
    task = swash.run(input_dir=file_manager.get_root_dir(),
                    sim_config_filename="input.sws",
                    on=machine_group)
```

The template mechanism will allow you to explore 5 different variations of the
simulation, each with a different water level. The simulations will be submitted
to our dedicated machine group and will run in parallel.
We can check that all simulations are running via the CLI and that it took only
1min for the moment they are submitted until they start running:

```bash
$ inductiva tasks list

    ID                         Simulator    Status    Submitted         Started           Computation Time    Resource Type

    57mr4kas99jxb9titkeackano  swash        started   01 Feb, 09:07:19  01 Feb, 09:08:03  *0:03:12            c2-standard-30
    ox8718m0pwfi02zczui3qky4w  swash        started   01 Feb, 09:07:17  01 Feb, 09:08:02  *0:03:14            c2-standard-30
    mak1ji62s7axf7mespkc36g7e  swash        started   01 Feb, 09:07:15  01 Feb, 09:08:03  *0:03:14            c2-standard-30
    ijyu8bkvme7vg9k0kj6v23gxa  swash        started   01 Feb, 09:07:14  01 Feb, 09:08:02  *0:03:16            c2-standard-30
    g5qq5c9mk2nr5wqhzef38sdm4  swash        started   01 Feb, 09:07:12  01 Feb, 009:08:01  *0:03:17            c2-standard-30
```

This is a great way to speed up the execution of multiple simulations, since the
time to run all 5 simulations will be approximately the same as running just one,
that is the above 5 simulations took 9m55s to complete, which is the time of
the slowest simulation.

Now, that all the simulations have finished running, we end this tutorial with an
extra lesson to help reduce the amount of time machines are left unused:

> Don't forget to terminate your computational resources with `inductiva resources terminate --all`.
