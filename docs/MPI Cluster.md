# Setting up a dedicated MPI Cluster on the cloud
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

## `MPICluster` Class

To instantiate an `MPICluster` object the following parameters can be configured:
- the `machine_type` defines the type of CPU used for each machine. This parameter
follows the naming convention set by [Google Cloud](https://cloud.google.com/compute/docs/machine-types),
e.g., `c2-standard-16`. This convention is composed of a prefix that defines the
CPU series, a suffix that sets the number of [virtual CPUs (vCPU)](https://cloud.google.com/compute/docs/cpu-platforms)
per machine and the middle word refers to the level of RAM per vCPU. In the example,
`c2` refers to an Intel Xeon Scalable processor of 2nd generation, `standard`
means 4 GB of RAM per vCPU and will contain `16` vCPUs.
Currently, this is the [list of available machine types available via the API]().
- the `num_machines` sets the number of machines available in the cluster. While the computational resource is active, these machines will be reserved
for the user.
- the `data_size_gb` allows the selection of the size of the disk attached to
each machine that is reserved for the simulation data in GB.

For example, the following code creates an MPICluster with 2 machines of type
`c2-standard-30`:

```python
import inductiva

mpi_cluster = inductiva.resources.MPICluster(
   machine_type="c2-standard-30",
   num_machines=2,
   data_size_gb=100)
```

When initializing an MPI cluster he is registered on the API, but the computational
resources won't be active yet. These can be launched with `mpi_cluster.start()`.
Within a few minutes, the machines will be set up and ready to collaborate on running simulations. At any moment, you can check an estimate of the price per
hour of the cluster with `mpi_cluster.estimate_cloud_cost()` and when you have finished
you can terminate it with `mpi_cluster.terminate()`. Running simulations will be killed and from this point, the `mpi_cluster` object cannot be re-used.

However, as you have seen, it is simple and clear how to instantiate a dedicated
MPI cluster on-demand without much hassle.

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