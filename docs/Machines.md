# Manage Computational Resources

By default, simulation requests are processed by a shared pool of machines serving 
multiple users. These machines are equipped with 4 cores each, and the pool is 
limited to a predefined capacity. If you require dedicated resources with custom 
configurations, the **Inductiva API** allows you to easily set up virtual machines 
reserved solely for your simulations. These are managed via machine groups, *i.e.*, 
groups of homogeneous machines with specific requested properties that can be 
started and terminated on demand. With the `MachineGroup` class, you can configure, 
start, and terminate machines. Creating a `MachineGroup` requires specifying the 
type of Virtual Machine to use and the number of such VMs that will compose the 
group. Currently, the available options for the `machine_type` are the ones 
available in the [Google Cloud Platform](https://cloud.google.com/compute/docs/machine-types). 
Once a `MachineGroup` is created, simply pass it as argument to your simulations, 
which will then be scheduled to run on those machines. 

Note that a `MachineGroup` is literally a group of individual machines that do 
not communicate with each other. In other words, a `MachineGroup` is not a 
computational cluster where the load of each simulation is divided over all 
machines of the cluster.

### Examples


Prepare a simulation, for example with Reef3D, and run it in a dedicated machine group:

```python

import inductiva

# Create a MachineGroup object with a single machine of the desired type
mg = inductiva.resources.MachineGroup(
    machine_type="c2-standard-16",
    num_machines=1,
    disk_size_gb=60,
)

# Start the machine group
mg.start()

# Download the configuration files for Reef3D
input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "reef3d-input-example.zip", unzip=True)

# Initialize the Reef3D simulator
simulator = inductiva.simulators.REEF3D()

# Run the simulation and wait for it to finish
task = simulator.run(input_dir=input_dir, on=mg)
task.wait()
```

Don't forget to always terminate your machine groups when you don't need them
anymore. If you want it to be terminated automatically after the simulation
completes, add at the end of the above code snippet:

```python
mg.terminate()
```

Otherwise, you can use the CLI to terminate it whenever you want:

```bash
inductiva machines terminate <machine_group_name>
```

#### Launch many simulations to run in parallel in a group of machines


```python

import inductiva

# Create and start a MachineGroup with 5 machines
mg = inductiva.resources.MachineGroup(
    machine_type="c2d-standard-16",
    num_machines=5)
mg.start()

# Download the configuration files for Reef3D
input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "reef3d-input-example.zip", unzip=True)

# Initialize the Reef3D simulator
simulator = inductiva.simulators.REEF3D()

tasks = []

# Submit 5 simulations at the same time to the same machine group.
# Each simulation will be picked up by a different machine of the group and all
# will run in parallel.
for _ in range(5):
    task = simulator.run(input_dir=input_dir, on=mg)
    tasks.append(task)

# Block until all tasks complete
for task in tasks:
    task.wait()

# Once you don't need them anymore, terminate the machines
mg.terminate()
```

#### List active machines

The machines you have running can be listed via the CLI:

```bash
inductiva machines list
```

or via Python

```python
import inductiva

inductiva.resources.machine_groups.list()
```

and provide the following information:

```bash
Name       Machine Type    Elastic       Type   # machines    Disk Size in GB       Spot   Started at (UTC)
api-369qthz2285cdg6f9exgyom31    c2d-highmem-112      False        mpi            2                200      False   22 Jan, 16:39:12
api-6g1eix70k2ixxsjtroyhuymgp    c2d-highmem-112      False   standard            1                200      False   22 Jan, 16:38:45
api-m141gzk927xd74szrw3mwim5f      c2-standard-4      False   standard            1                 60      False   23 Jan, 10:20:59
api-mdam59yq17m9vsaktk8b98bdr      c2-standard-4      False   standard            1                 60      False   23 Jan, 10:22:43
```

Moreover, the active machines can be retrieved and used via Python, as follows:
```python
# Get a list of all the MachineGroup objects (for example, if you want to terminate them all at once)
mg_list = inductiva.resources.machine_groups.get()
mg_list
#[<inductiva.resources.machines.MachineGroup at 0x7f8cde53d2a0>,
# <inductiva.resources.machines.MachineGroup at 0x7f8c58954c70>]

# Terminate the first machine of the list
mg_list[0].terminate()
```
