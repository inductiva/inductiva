## Set up an Elastic Machine Group

The `ElasticMachineGroup`, similarly to the standard `MachineGroup`, is composed of 
a group of homogeneous machines that work individually to run multiple simulations. 
The difference is that the number of active machines is scaled up and/or down 
automatically based on the simulations in the queue. This prevents computational 
resources from being idle when there are no sufficient simulations to run and 
allows scaling up the computational resources when the queue is full.

Note that, the elasticity is independent of each machine being preemptible, i.e.,
these can be combined and Inductiva API manages the simulations accordingly.

Based on the parameters above you can configure it as you wish, and then start it 
to make it available for your simulations.

Let's start an `ElasticMachineGroup`:

```python
import inductiva

# Configure an elastic machine group to start with a minimum of 1 machine up to a
# maximum of 5, each with a disk
# size of 70 GB and an Intel Xeon Scalable processor of 2nd gen with 4 vCPUs and 
# 4 GB of RAM per vCPU.
elastic_machine_group = inductiva.resources.ElasticMachineGroup(
    machine_type="c2-standard-4",
    min_machines=1,
    max_machines=5,
    data_disk_gb=70,
    spot=False
)

# Launch the Elastic machine group to make it available to run simulations:
elastic_machine_group.start()
```

Once started, it can be passed as an argument on the simulator `run` method and it 
will be added to the queue of the computational resource.
