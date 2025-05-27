# ElasticMachineGroup Class

An Elastic machine group is similar to the [Machine group](./machinegroup_class.md)
with the extra property that scales up and down the number of active machines based
on the number of simulations in queue. It is composed of a pool of homogeneous machines
that work individually and do not communicate with each other in any way. 

Hence, an elastic machine group creates a private queue for which workers scale
based on the number of tasks in it. This allows running multiple simulations at the
same time, with the slight overhead of machines starting, with a more 
cost-effective strategy since machines won't stay idle for long.

## Instantiating an `ElasticMachineGroup` object
To create an elastic machine group the following properties can be configured:
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
[list of available machine types available via the API.](https://tutorials.staging.inductiva.ai/intro_to_api/computational-infrastructure.html#available-computational-resources)
- the `min_machines`, `max_machines` sets the number of minimum and maximum machines 
available in the computational resource. That is, the number of active machines will
never go lower than the minimum and never above the maximum. During runtime, there
might be a different number of active machines in between. Moreover, the `min_machines``
is the number of machines that the group is started.
- the `data_disk_gb` allows the selection of the size of the disk attached to each
machine that is reserved for the simulation data in GB.
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

elastic_machine_group = inductiva.resources.ElasticMachineGroup(
    machine_type="c2-standard-16",
    min_machines=2,
    max_machines=10,
    data_disk_gb=100,
    spot=False)
```

Creating an instance of `ElasticMachineGroup` does not start the machines. This only 
registers the configuration on the API which can now be used to manage it further.

## Managing the ElasticMachineGroup

With your `elastic_machine_group` object ready, you can launch the elastic machine
group with the minimum number of machines active with `elastic_machine_group.start()`.

Within a few minutes, the machines will be set up and ready to pick up several
simulations simultaneously. As simulations get into the queue, the number of active
machines increases.

At any moment, you can check an estimate of the price per hour of the group as follows:

```
elastic_machine_group.estimate_cloud_cost()
```

When you have finished you can terminate it with `elastic_machine_group.terminate()`
or via the CLI with `$ inductiva resources terminate api-agn23rtnv0qnfn03nv93nc`.

Running simulations will be killed and from this point, the `elastic_machine_group`
object cannot be re-used.

Elastic Machine Group on demand without any hassle.
