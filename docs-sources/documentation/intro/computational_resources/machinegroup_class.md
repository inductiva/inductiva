# MachineGroup Class

A Machine group is a pool of homogeneous machines that work individually and 
which do not communicate with each other in any way. Launching a machine group
allows the creation of a private queue that only receives the tasks you specifically
send to them. Then, the machines can pick simulations from the queue, which allows
to run multiple simulations in parallel and speeds up the exploration of a design space.

To instantiate a `MachineGroup` object the following parameters can be configured:
- the `machine_type` defines the type of CPU used for each machine. This parameter
follows the naming convention set by
[Google Cloud](https://cloud.google.com/compute/docs/machine-types),
e.g., `c2-standard-16`. This convention is composed of a prefix that defines the
CPU series, a suffix that sets the number of
[virtual CPUs (vCPU)](https://cloud.google.com/compute/docs/cpu-platforms)
per machine and the middle word refers to the level of RAM per vCPU. In the example,
`c2` refers to an Intel Xeon Scalable processor of 2nd generation, `standard`
means 4 GB of RAM per vCPU and will contain `16` vCPUs.
Currently, this is the 
[list of available machine types available via the API](https://tutorials.staging.inductiva.ai/intro_to_api/computational-infrastructure.html#available-computational-resources). 
- the `num_machines` sets the number of machines available in the computational
resource. While the computational resource is active, these machines will be reserved
for the user.
- the `data_disk_gb` allows the selection of the size of the disk attached to each machine that is reserved for the simulation data in GB.
- the `spot` argument determines if the machines will be preemptible or standard.
Preemptible machines can be stopped at any time and for that reason are only
advised for fault-tolerant workloads. If simulations are running when they are
stopped, the simulation is resubmitted to the queue of the machine group again.
- the `max_idle_time` determines the time a machine group can remain idle (without
receiving any task) before it is terminated.
- the `auto_terminate_ts` defines the moment in time in which the resource will
be automatically terminated, even if there are tasks still running.

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

Creating an instance of `MachineGroup` does not start the machines. 
This only registers the configuration on the API, which can now be used
to manage it further.

## Managing the MachineGroup

With your `machine_group` object ready, starting all of the machines at the same
time is as simple as calling `machine_group.start()`.

Within a few minutes, the machines will be set up and ready to pick several
simulations simultaneously. At any moment, you can check an estimate of the
price per hour of the group with `machine_group.estimate_cloud_cost()` and
when you have finished you can terminate it with `machine_group.terminate()`.
Running simulations will be killed and from this point, the `machine_group`
object cannot be re-used.

To simplify the workflow, the last two functions can also be performed via the CLI.

First, you can check the cost of the group by selecting the machine
type and the number of machines you wish to use:

```bash
$ inductiva resources cost c2-standard-4 -n 4
Estimated total cost (per machine): 0.919 (0.230) $/h.
```

When you don't need the Machine group anymore, you can easily kill it with the name:

```bash
$ inductiva resources terminate api-agn23rtnv0qnfn03nv93nc
```

Machine Group on demand without any hassle.
