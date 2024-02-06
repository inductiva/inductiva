# Overview

In this guide, we will explain some of the main features of the Inductiva API when 
it comes to making informed decisions about resource allocation. Here, you will 
learn the two options available for running your simulations; shared and dedicated
resources through the Inductiva API. Lastly, we 
demonstrate the performance differences between shared and dedicated resources 
using a SWASH simulation.

## Shared Resources

As a standard practice, when you submit simulation tasks through the Inductiva API, 
they enter a [task queue]() and are sent to a shared pool of workers on designated Virtual Machines (VMs)
serving multiple users from either the [Google Cloud Provider (GCP)](https://cloud.google.com/compute/docs/machine-resource) or Inductiva's own 
computational platform (ICE). These VMs are specifically allocated to facilitate easy 
API testing and to handle light tasks efficiently with a minimal setup, ideal for 
quick experimentation.

However, this shared resource pool has its limitations, including slower completion 
times for simulations due to its finite capacity and less powerful VMs. 

For running a larger volume of simulations that require more powerful VMs, the 
Inductiva API provides you with the option to set up your own [dedicated resources](). 
to ensure you can access the necessary computational power for your simulations 
without the constraints of a shared environment. 

## Dedicated Resources

Inductiva provides the capability to create dedicated pools of VM resources, 
termed [Machine Groups](), exclusively reserved for your use and not shared with 
others. 

In this option, users can launch three types of computational resources for their simulations:

- [**Machine Group**](#launch-a-machine-group): This consists of homogeneous machines 
designed to operate individually, enabling the distribution of multiple simulations 
across different machines for parallel processing.
- [**Elastic Machine Group**](#set-up-an-elastic-machine-group): Similar to Machine 
Group, these also consist of individual machines. The key advantage here is the 
elastic scaling feature, which dynamically adjusts the number of machines based 
on simulation demands, ensuring efficient resource utilization.
- [**MPI Cluster**](#start-a-mpi-cluster-in-the-cloud) This setup involves a network 
of machines configured to work in tandem on a single simulation task, distributing 
the workload across multiple CPUs. This is particularly useful for complex simulations 
that exceed the capabilities of a single machine.


## Comparison Example: SWASH Simulation

In this comparison, we aim to illustrate the performance differences between 
running simulations on a shared pool of resources and on dedicated resources using the 
Inductiva API. 

We will use the [SWASH simulation]() as our example to demonstrate these differences.

### SWASH Simulation on Shared Resources

First, we'll run the simulation using the shared pool of workers, a convenient 
option for those getting started or running less resource-intensive tasks. 

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
In this option, it takes **approximately 20 minutes** 
to run this simulation on the shared pool of resources.

### SWASH Simulation on Dedicated Resources

Let's now run the same simulation on dedicated resources, specifically set 
up for this task.

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
Notice that, the simulation is picked almost immediately - no waiting time required, and by contrasting the two resource options above, running the simulation on a [dedicated machine group]() with a `c2-standard-30` machine took **9 minutes and 37s** to complete, which is 
**2.68 times less** than on the [shared pool](). 

## What to read next

Learn how to [run multiple simulations in parallel]() 

Learn how to [set up an API cluster]()









