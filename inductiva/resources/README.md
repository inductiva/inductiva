# Manage Computational Resources

By default, simulation requests are processed by a shared pool of machines serving multiple users. These machines are equipped with 4 cores each, and the pool is limited to a predefined capacity. If you require dedicated resources with custom configurations, the **Inductiva API** allows you to easily set up virtual machines reserved solely for your simulations. These are managed via machine groups, *i.e.*, groups of homogeneous machines with specific requested properties that can be started and terminated on demand.

There are two types of machine groups available: `MachineGroup` and `ElasticMachineGroup`. Both of these classes serve to configure, start, and terminate the machines, and will require specifying the type of Virtual Machine to use. At the moment only [general-purpose](https://cloud.google.com/compute/docs/general-purpose-machines), [compute optimized](https://cloud.google.com/compute/docs/compute-optimized-machines), and [memory optimized](https://cloud.google.com/compute/docs/memory-optimized-machines) Google Cloud `machine_types` are available. The only difference is in a way these groups manage the number of machines up.

- `MachineGroup` creates the specified number of machines at once, those machines will be up until they are terminated.
- `ElasticMachineGroup` starts with a minimum number of machines and dynamically scales up to a maximum number of machines or down to minimum based on the CPU load.

Once a machine group is created, simply pass it as argument to your simulations, which will then be scheduled to run on those machines.

Note that these machine groups are literally groups of individual machines that do not communicate with each other. In other words, a machine group is not a computational cluster where the load of each simulation is divided over all machines of the cluster.

### Examples


#### Running a simulation in a specific machine type
#### Launch `MachineGroup`

```python

import inductiva
from inductiva import molecules

# Create a MachineGroup object with a single machine of the desired type
mg = inductiva.resources.MachineGroup(
    machine_type="c2-standard-4",
    num_machines=1,
    disk_size_gb=60,
)

# Estimate the hourly cost (in $) of the cloud resources requested
price_per_hour = mg.estimate_cloud_cost()

# Start the machines
mg.start()
```
#### Launch `ElasticMachineGroup`

```python

import inductiva
from inductiva import molecules

# Create a ElasticMachineGroup object with minimum 1 machine. Based on the CPU load this machine group
# will scale up to 5 machines, if there is no load it will scale down to 1 machine.
mg = inductiva.resources.ElasticMachineGroup(
    machine_type="c2-standard-4",
    min_machines=1,
    max_machines=5,
    disk_size_gb=60,
)

# Estimate the hourly cost (in $) of the cloud resources requested
price_per_hour = mg.estimate_cloud_cost()

# Start the machines
mg.start()
```

Create your simulation scenario and run it on your machine group. Example with ProteinSolvation scenario:

```python

# Download the insulin protein (ID - "1ZNI") from RCSB database
insulin_pdb_file = inductiva.molecules.utils.download_pdb_from_rcsb(pdb_id="1ZNI")

scenario = molecules.ProteinSolvation(
    protein_pdb=insulin_pdb_file, temperature=300)

# Pass your machine group object when submitting a simulation so that it runs
# on the desired machine group
task = scenario.simulate(machine_group=mg, ignore_warnings=True)


# Once your simulations are done, terminate the machines
mg.terminate()
```

Once the machine group with the desired configurations are started, you can run hundreds of simulations in parallel in those machines.


#### List and get active machine groups

You can also list your active machine groups, or get a list of `MachineGroup` objects of previously created machine groups:

```python
import inductiva

inductiva.resources.machine_groups.list()
# INFO:absl:Active machine groups:
#                                     Name         VM Type   # machines    Disk Size in GB       Spot         Started at
# api-1b1f724c-5cfe-4d87-8439-9689aa139723   c2-standard-4          1/1                 60      False   13 Sep, 07:38:50
# api-8e6bf7d8-4888-4de9-bda5-268484b46e6f   c2-standard-4          2/5                 60      False   13 Sep, 07:37:49


# Get a list of all the MachineGroup objects (for example, if you want to terminate them all at once)
mg_list = inductiva.resources.machine_groups.get()
mg_list
#[<inductiva.resources.machines.MachineGroup at 0x7f8cde53d2a0>,
# <inductiva.resources.machines.MachineGroup at 0x7f8c58954c70>]
```


## Resources Quotas

As a user, you are allowed to use resources up to certain limits for free. The current free usage limits are:

- Up to 80 cores in total at a time, i.e. the sum of all active cores should not exceed 80
- Up to 10 machines at a time (example: you can have 2 machine groups with 5 machines each, or 1 machine group with 10 machines)
