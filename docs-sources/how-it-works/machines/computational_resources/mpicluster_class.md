# MPICluster Class

An MPI cluster is a set of machines that work together to solve a problem. These
machines are configured to communicate with each other using the [Message Passing
Interface (MPI)](https://en.wikipedia.org/wiki/Message_Passing_Interface) protocol.

With the `MPICluster` class, you can create a dedicated MPI cluster on-demand that
is configured automatically to run your simulations in a distributed manner across
the machines.

## Instantiating an `MPICluster` object
The following parameters can be configured:
- the `machine_type` defines the type of CPU used for each machine. This parameter
follows the naming convention set by
[Google Cloud](https://cloud.google.com/compute/docs/machine-types),
e.g., `c2-standard-16`. This convention is composed of a prefix that defines the
CPU series, a suffix that sets the number of
[virtual CPUs (vCPU)](https://cloud.google.com/compute/docs/cpu-platforms)
per machine and the middle word refers to the level of RAM per vCPU. In the example,
`c2` refers to an Intel Xeon Scalable processor of 2nd generation, `standard`
means 4 GB of RAM per vCPU and will contain `16` vCPUs.
Check out the 
[complete machine catalog available via the API](https://inductiva.ai/machines).
- the `zone` allows to select the zone where the machines will be launched. By default, machines are launched in the `europe-west1-b` zone.
- the `num_machines` sets the number of machines available in the cluster.
While the computational resource is active, these machines will be reserved
for the user.
- the `data_disk_gb` defines the size of the NFS partition that is mounted on
the head node and shared with the worker nodes. This partition is used to store
the input and output files of the simulations.
- the `max_idle_time` determines the time a machine group can remain idle (without
receiving any task) before it is terminated. By default, this value is _3 minutes_.
- the `auto_terminate_ts` defines the moment in time in which the resource will
be automatically terminated, even if there are tasks still running.

For example, the following code creates an MPICluster with 2 machines of type
`c2-standard-30`:

```python
import inductiva

mpi_cluster = inductiva.resources.MPICluster(
   machine_type="c2-standard-30",
   num_machines=2,
   data_disk_gb=100)
```

Creating an instance of `MPICluster` does not start the machines. This only registers
the configuration on the API which can now be used to manage the cluster further.

## Managing the MPI Cluster

Visit our [Manage Resources](../manage_computational_resources.md) guide to learn how to monitor and control your `MPICluster` resources.
