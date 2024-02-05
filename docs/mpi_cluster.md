# Set up a dedicated MPI Cluster
In this guide, we extend the computational resources to empower the user ability
to run large scaling simulations on a dedicated MPI Cluster.

### What We'll Cover
* [Overview]()
* [Why will you need an MPICluster?]()
* [MPICluster Class]()
* [Running a simulati]
* [What to read next]()

## Overview
>*we need something here, I'll think it through - Maya*

## Why will you need an MPICluster?

At times, when simulations grow in complexity and size they may not fit
on a single machine, whether because the RAM is not enough or the simulation takes
too much time to finish running over a small number of cores.

In those cases, the simulations will need to be run on a cluster of machines that
work together to solve the problem. This is where the MPICluster comes in.

[Message Passing Interface](https://en.wikipedia.org/wiki/Message_Passing_Interface),
more well-known as MPI, is a standard for message-passing data on parallel computing architectures, in particular, over multiple machines. It is widely used in
the simulation world, and most open-source simulators integrated within the Inductiva API use it to scale their simulations.

However, setting up the machines to communicate with each other in a way that MPI
can thereafter take the stage, is no easy task and requires expertise to obtain
a coherent and efficient setup.

With Inductiva API you won't need to be an expert in setting up an MPI Cluster, we
have already taken care of that for you. In the following, you will learn how
to launch an MPI cluster to the cloud with only a few lines of code and run your
simulations over multiple machines.

Recall, that the SWASH simulation we ran in a [previous tutorial](Machines.md)
took 25m37s on the shared pool of workers and 9m37s on a selected `MachineGroup`
with a `c2-standard-30` machine type. We will use this as a baseline to compare
with the performance of the MPI Cluster.
 
> TODO LUIS

### Distributing a simulation over multiple machines

#### Example

Now that we explained what this class is about, let’s see how to use it to run
a SWASH simulation and understand the gains of using an `MPICluster`.

```python
import inductiva

# Download the input files for the SWASH simulation
input_dir = inductiva.utils.download_from_url(
   "https://storage.googleapis.com/inductiva-api-demo-files/"
   "swash-resources-example.zip", unzip=True)

# Instantiate a MPICluster object with 4 machine of type c2-standard-30 and start it
# immediately. This accounts for 120 vCPUs
mpi_cluster = inductiva.resources.MPICluster(
   machine_type="c2-standard-30", num_machines=4)
mpi_cluster.start()

# Initialize the SWASH simulator and run the simulation
# in your just launched dedicated MPICluster
swash = inductiva.simulators.SWASH()
task = swash.run(input_dir=input_dir,
                sim_config_filename="input.sws",
                on=mpi_cluster)

# Wait for the task to finish and download the outputs
task.wait()

# Terminate your dedicated MPICluster at then end of the simulation.
mpi_cluster.terminate()
```

With the MPI cluster using 120 vCPUs, the simulation took 3m25s to finish, which
is a 2.75 times reduction compared to using a single machine with 30 vCPUs.
Notice that the time reduction is not linear with the number of vCPUs, but it
is still a significant improvement. These speed-up are significant when running
longer simulations and using a readily available MPI cluster can reduce the
time to obtain results from one day to a few hours.


## What to read next
* [Set up an MPI]()
* [Get an overview of the CLI]()