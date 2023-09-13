# Manage Computational Resources

By default, simulation requests are processed by a shared pool of machines serving multiple users. If you require dedicated resources, the **Inductiva API** allows you to easily set up virtual machines reserved solely for your simulations. These are managed via machine groups, *i.e.*, groups of homogeneous machines with specific requested properties that can be started and terminated on demand. With the `MachineGroup` class, users can configure, start, and terminate machines. Creating a `MachineGroup` requires specifying the type of Virtual Machine to use and the number of such VMs that will compose the group. Currently, the available options for the `machine_type` are the ones available in the [Google Cloud Platform](https://cloud.google.com/compute/docs/machine-types). Once a `MachineGroup` is created, simply pass it as argument to your simulations, which will then be scheduled to run on those machines. Note that a `MachineGroup` is literally a group of individual machines that do not communicate with each other. In other words, a `MachineGroup` is not a computational cluster where the load of each simulatin is divided over all machines of the cluster

### Examples


#### Running one simulation in a specific machine type

```python

import inductiva
from inductiva import molecules

# Create a MachineGroup object with a single machine of the desired type
mg = inductiva.resources.MachineGroup(
    machine_type="c2-standard-4",
    num_machines=1,
    disk_size_gb=40,
)

# Estimate the hourly cost (in $) of the cloud resources requested
price_per_hour = mg.estimate_cloud_cost()

# Start the machines
mg.start()

# Example with ProteinSolvation scenario
scenario = molecules.ProteinSolvation(
    protein_pdb=my_pdb_file, temperature=300)

# Pass your machine group object when submitting a simulation so that it runs
# on the desired machine
task = scenario.simulate(machine_group=mg, run_async=False)


# Once you don't need it anymore, terminate the machines
mg.terminate()
```

#### Launch hundreds of simulations to run in parallel in a group of machines


```python

import inductiva
from inductiva import molecules

# Create a MachineGroup with 10 machines
mg = inductiva.resources.MachineGroup(
    machine_type="c2d-standard-16",
    num_machines=10,
    disk_size_gb=40,
)

price_per_hour = mg.estimate_cloud_cost()

# Start the machines
mg.start()

# Create a list of e.g. 100 protein files you want to simulate
my_pdb_files = get_my_pdb_file_paths()
tasks = []

# And simulate the solvation of the proteins on your newly started machines
for my_pdb_file in my_pdb_files:
    scenario = molecules.ProteinSolvation(
        protein_pdb=my_pdb_file, temperature=300)

    # Tasks will be submitted asynchronously to the new machine group
    task = scenario.simulate(machine_group=mg, run_async=True)
    tasks.append(task)

# Block until all tasks complete
for task in tasks:
    task.wait()

# Once you don't need them anymore, terminate the machines
mg.terminate()
```

#### List active machine groups

You can also list your active machine groups, and use the names that appear in the console to recreate a `MachineGroup` object of a previously created machine group:

```python
import inductiva

inductiva.resources.list_active_machine_groups()
#                                     Name         VM Type   # machines    Disk Size in GB       Spot         Created at
# api-1b1f724c-5cfe-4d87-8439-9689aa139723   c2-standard-4            1                 40      False   13 Sep, 07:38:50
# api-8e6bf7d8-4888-4de9-bda5-268484b46e6f   c2-standard-4            1                 40      False   13 Sep, 07:37:49

# Create a MachineGroup object by using its name and resume it using
mg = inductiva.resources.get_machine_group("api-24e497af-d135-4a59-bdd1-854bf0176cbf")

# Or get a list of all the MachineGroup objects (for example, if you want to terminate it at once)
mg_list = inductiva.resources.get_machine_groups()
```