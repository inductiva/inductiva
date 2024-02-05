# `MachineGroup` Class

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