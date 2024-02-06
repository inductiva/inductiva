
> TO BE BROKEN INTO PIECES

# Computational Resources

Besides running and managing simulations, Inductiva API enhances the computational 
infrastructure of any user by providing a simple interface to run simulations in
your own dedicated computational resources.

In this tutorial, we learn the API capabilities to launch and manage 
computational resources in the cloud and guide the scaling of your 
computational infrastructure to run simulations faster and more efficiently.

Before getting into the details, we highlight that launching a computational 
resource and running a simulation in it is as simple as follows:

```python
import inductiva

# Configure a computational resources
machines = inductiva.resources.MachineGroup(
    machine_type="c2-standard-16", num_machines=5, data_disk_gb=60)

# Launch the computational resource to be available to run simulations
machines.start()

# Run the Getting Started example simulation in our own dedicated computational
# resource
input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "reef3d-dambreak-dir.zip", unzip=True)

simulator = inductiva.simulators.REEF3D()

# Select the computational resource to run the simulation on the `run` method
task = simulator.run(input_dir=input_dir, on=machines)

task.wait()

# Terminate the computational resource when you don't need it anymore
machines.terminate()
```

This code snippet launches a computational resource with 5 machines of type
`c2-standard-16` that are available to now deploy your own simulations. In the
following sections, we will dive deep into the details of the computational resources
available via the API and how to use them.

## Running simulations in parallel

#### Example

The second use case of launching a `MachineGroup` is that of setting multiple
simulations running in parallel and distributed by the various machines constituting
the group. This is useful when you want to run multiple simulations in parallel,
but you don't want to wait for the first one to finish before starting the second one. 

To exemplify, we will use the [templating mechanism]() built-in the Inductiva API
to automatically change the water level of the simulation in the input files and
run 5 different simulations in parallel. 

```python
import inductiva
from inductiva import mixins

# Instantiate a MachineGroup object with 1 preemptible machine of type
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
-------------------------  -----------  --------  ----------------  ----------------  ------------------  ---------------
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

#### Configuration parameters

The above computational resources are configured with a few common parameters. Let's
start by introducing them:

- `machine_type` defines the type of CPU used for each machine. This parameter
follows the naming convention set by [Google Cloud](https://cloud.google.com/compute/docs/machine-types).
Each machine type, e.g. `c2-standard-60`, is composed of a prefix defining the
CPU series, a suffix that sets the number of [virtual CPUs (vCPU)](https://cloud.google.com/compute/docs/cpu-platforms)
per machine and the middle word refers to the level of RAM per vCPU. In the example,
`c2` refers to an Intel Xeon Scalable processor of 2nd generation, `standard`
means 4 GB of RAM per vCPU and will contain `60` vCPUs. Below, we introduce the
machine types available via the API. 
- the parameters `num_machines`, `min_machines`, `max_machines` control the number
of machines available in the computational resource. The former is used by the
Machine group and the MPI Cluster resources. The latter ones set the minimum number
of machines that are always available and the maximum number of machines that can
be available at any time, respectively. 
- `data_disk_gb` allows the selection of the size of the disk attached to each machine that is reserved for the simulation data in GB.
- `spot` allows to differ between standard resources or preemptible ones. The latter
resources are way cheaper but only serve for fault-tolerant workloads since they
can be stopped at any time. The former are fully dedicated to the user's usage.
In case of failure of a simulation with a preemptible machine, the simulation is
resubmitted to the queue of the computational resource. Currently, this is only
available for the Machine Group and Elastic Machine Group.



