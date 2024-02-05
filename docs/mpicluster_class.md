# `MPICluster` Class

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
- the `data_disk_gb` allows the selection of the size of the disk attached to
each machine that is reserved for the simulation data in GB.

For example, the following code creates an MPICluster with 2 machines of type
`c2-standard-30`:

```python
import inductiva

mpi_cluster = inductiva.resources.MPICluster(
   machine_type="c2-standard-30",
   num_machines=2,
   data_disk_gb=100)
```

When initializing an MPI cluster he is registered on the API, but the computational
resources won't be active yet. These can be launched with `mpi_cluster.start()`.
Within a few minutes, the machines will be set up and ready to collaborate on running simulations. At any moment, you can check an estimate of the price per
hour of the cluster with `mpi_cluster.estimate_cloud_cost()` and when you have finished
you can terminate it with `mpi_cluster.terminate()`. Running simulations will be killed and from this point, the `mpi_cluster` object cannot be re-used.

However, as you have seen, it is simple and clear how to instantiate a dedicated
MPI cluster on-demand without much hassle.