# MachineGroup Class

A Machine Group is a pool of homogeneous machines that work individually and 
which do not communicate with each other in any way. Launching a Machine Group
allows the creation of a private queue that only receives the tasks you specifically
send to them. Then, the machines can pick simulations from the queue, which allows
to run multiple simulations in parallel and speeds up the exploration of a design space.

## Instantiating a `MachineGroup` object
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
- the `num_machines` sets the number of machines available in the computational
resource. While the computational resource is active, these machines will be reserved
for the user.
- the `data_disk_gb` allows the selection of the size of the disk attached to each machine that is reserved for the simulation data in GB.
- the `spot` argument determines if the machines will be preemptible or standard.
Preemptible machines can be stopped at any time and for that reason are only
advised for fault-tolerant workloads. If simulations are running when they are
stopped, the simulation is resubmitted to the queue of the machine group again.
- the `max_idle_time` determines the time a machine group can remain idle (without
receiving any task) before it is terminated. By default, this value is _3 minutes_.
- the `auto_terminate_ts` defines the moment in time in which the resource will
be automatically terminated, even if there are tasks still running.

For example, the following code creates a `MachineGroup` with 2 machines of type
`c2-standard-16` with 100 GB of disk space each:

```python
import inductiva

machine_group = inductiva.resources.MachineGroup(
    machine_type="c2-standard-16",
    num_machines=2,
    data_disk_gb=100,
    spot=False)

machine_group.start()  # start the MachineGroup
```

Creating an instance of `MachineGroup` does not start the machines. 
This only registers the configuration on the API, which can now be used
to manage it further.

## Managing the Machine Group

Visit our [Manage Resources](../manage_computational_resources.md) guide to learn how to monitor and control your `MachineGroup` resources.

```{banner_small}
:origin: how_it_works_mg_class
```