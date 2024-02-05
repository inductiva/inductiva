# Simulation Processing and Computational Resource Allocation
In this guide, we will explain some of the main features of the Inductiva API when 
it comes to managing and utilizing computational resources.

### What We'll Cover
* [Overview]()
* [How can you set up computational resources with Inductiva?]()
    * [Task Queue and Shared VM Pool Integration]()
    * [Custom Hardware Setup for Enhanced Simulation Performance]()
* [MachineGroup Class]()
    * Setting up a MachineGroup for selecting specific VM type
    * Setting up a MachineGroup for running simulations in parallel
* [What to read next]()


## Overview
>*we need something here, I'll think it through - Maya*

## How can you set up computational resources with Inductiva?

>*an intro to the two options - Maya*
### Task Queue and Shared VM Pool Integration

Simulation tasks you create and submit via the Inductiva will be sent to a Task 
queue, and will eventually be picked up by a pool of workers running on Virtual 
Machines (VMs) available from a Cloud provider

By default, simulation tasks will be sent to a shared pool of workers serving 
multiple users. These workers live on VMs that we decided to set aside to make it 
easier for any user to test the API, and run also relatively light tasks with 
multiple users. These workers live on VMs that we decided to set aside to make it 
easier for any user to test the API, and run also relatively light tasks with 
the simplest possible setup. 

For example, the code below will start a SWASH simulation that will be 
automatically picked up by the shared pool of workers.

**WARNING:** For the sake of demonstrating performance differences, the following
simulation takes around 20 minutes to complete.

```python
import inductiva

# Download the input files for the SWASH simulation
input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "swash-resources-example.zip", unzip=True)

# Initialize the SWASH simulator and run the simulation
# in the default shared pool of workers
swash = inductiva.simulators.SWASH()
task = swash.run(input_dir=input_dir,
                 sim_config_filename="input.sws")

# Wait until the task finishes running.
task.wait()
```

Observe that at no point we explicitly defined the target VMs where this simulation
would be executed, or even just their specs. Instead, the task will get automatically
sent to the shared pool of workers that we prepared for all users. 
This is very simple, and a great way for doing quick experimentation.  

However, despite the convenience and simplicity, the above simulation took 25m37s
to complete. The shared pool of resources 
has a limited predefined capacity and doesn't possess powerful VMs. Therefore, since
it is shared by all users, it is not appropriate for executing larger tasks,
since waiting times can be extremely large. So, if you need to run a larger number
of simulation tasks, and you need more powerful VMs to run it, you will need to
reserve that capacity for your exclusive use.

### Custom Hardware Setup for Enhanced Simulation Performance

Inductiva provides a way of creating pools of VMs resources that are exclusively 
available for you, and not shared with any other user. We call these Machine Groups, 
and, as we will explain later, they come in two flavors. 

A Machine Group are groups of homogeneous cloud VMs with specific specs that you 
can define programmatically and terminate on demand via the API. This will give 
you full control of the type of VM you use to run your simulations, and will ensure 
a certain amount of compute power that we reserve exclusively for you. 

Note that a Machine Group is literally a group of individual VMs that do not 
communicate with each other. In other words, a MachineGroup is not a cluster, 
such as an MPI Cluster, where the load of each simulation is divided over all 
machines of the cluster. To set up an MPI Cluster, see here (point to the MPI 
documents). 

But for now, let’s dig deeper in the MachineGroup class that implements this 
notion of a group of independent machines of the same type. 
 
## `MachineGroup` Class

To instantiate a `MachineGroup` object the following parameters can be configured:
- the `machine_type` defines the type of CPU used for each machine. This parameter
follows the naming convention set by [Google Cloud](https://cloud.google.com/compute/docs/machine-types),
e.g., `c2-standard-16`. This convention is composed of a prefix that defines the
CPU series, a suffix that sets the number of [virtual CPUs (vCPU)](https://cloud.google.com/compute/docs/cpu-platforms)
per machine and the middle word refers to the level of RAM per vCPU. In the example,
`c2` refers to an Intel Xeon Scalable processor of 2nd generation, `standard`
means 4 GB of RAM per vCPU and will contain `16` vCPUs.
Currently, this is the [list of available machine types available via the API]().
- the `num_machines` sets the number of machines available in the computational
resource. While the computational resource is active, these machines will be reserved
for the user.
- the `data_disk_gb` allows the selection of the size of the disk attached to each machine that is reserved for the simulation data in GB.
- the `spot` argument determines if the machines will be preemptible or standard.
Preemptible machines can be stopped at any time and for that reason are only
advised for fault-tolerant workloads. If simulations are running when they are
stopped, the simulation is resubmitted to the queue of the machine group again.

For example, the following code creates a MachineGroup with 2 machines of type
`c2-standard-16` with 100 GB of disk space each:

```python
import inductiva

machine_group = inductiva.resources.MachineGroup(
    machine_type="c2-standard-16",
    num_machines=2,
    data_disk_gb=100,
    spot=False)
```

The instantiation of one of these objects registers the configuration on the API,
but no resources are active yet. These can be launched with `machine_group.start()`.
Within a few minutes, the machines will be ready to use and thereafter you can launch
your simulations there. At any moment, you can check an estimate of the price per
hour of the machine group with `machine_group.estimate_cloud_cost()`.
When you are done with the machines, you can terminate them with `machine_group.terminate()`.
Running simulations will be killed. From this point, the `machine_group` cannot be
re-used. But as you have seen it is simple to just instantiate a new one.

Now that we explained what this class is about, let’s see how to use it for two 
different purposes. First, we will create a MachineGroup with 1 machine only. 
The goal is to use this functionality to select specific types of VM available on 
Google Cloud for running our simulation. We will show the impact on performance 
(and potential cost) of running the simulation on different types of VMs.

Then, we will create a MachineGroup with 5 instances and show how to run 5 
variations of the same simulation in parallel. 

### Setting up a Machine Group for selecting specific VM type

#### Example

Let's now run the above simulation in our own dedicated resource.

```python
import inductiva

# Download the input files for the SWASH simulation
input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "swash-resources-example.zip", unzip=True)

# Instantiate a MachineGroup object with 1 preemptible machine of type
# c2-standard-30 and start it immediately
machine_group = inductiva.resources.MachineGroup(
    machine_type="c2-standard-30", spot=True)
machine_group.start()

# Initialize the SWASH simulator and run the simulation
# in your just launched dedicated MachineGroup
swash = inductiva.simulators.SWASH()
task = swash.run(input_dir=input_dir,
                 sim_config_filename="input.sws",
                 on=machine_group)

# Wait for the task to finish and download the outputs
task.wait()

# Terminate your dedicated MachineGroup at then end of the simulation.
machine_group.terminate()
```

Running the same simulation on a dedicated machine group with a `c2-standard-30`
machine took 9m37s, which is 2.68 times less than on the shared pool. Notice
that, the simulation is picked almost immediately - no waiting time required - and
selecting a more powerful machine greatly reduced the execution time.

### Setting up a MachineGroup for running simulations in parallel

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


## What to read next
* [Set up an MPI]()
* [Get an overview of the CLI]()
