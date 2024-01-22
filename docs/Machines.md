# Manage Computational Resources

By default, simulation requests are processed by a shared pool of machines serving multiple users. These machines are equipped with 4 cores each, and the pool is limited to a predefined capacity. If you require dedicated resources with custom configurations, the **Inductiva API** allows you to easily set up virtual machines reserved solely for your simulations. These are managed via machine groups, *i.e.*, groups of homogeneous machines with specific requested properties that can be started and terminated on demand. With the `MachineGroup` class, you can configure, start, and terminate machines. Creating a `MachineGroup` requires specifying the type of Virtual Machine to use and the number of such VMs that will compose the group. Currently, the available options for the `machine_type` are the ones available in the [Google Cloud Platform](https://cloud.google.com/compute/docs/machine-types). Once a `MachineGroup` is created, simply pass it as argument to your simulations, which will then be scheduled to run on those machines. 

Note that a `MachineGroup` is literally a group of individual machines that do not communicate with each other. In other words, a `MachineGroup` is not a computational cluster where the load of each simulation is divided over all machines of the cluster.

## Examples


### Running a simulation in a specific machine type

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
```

Launch your simulation in the machine group you just launched. For example:
```python
# Set simulation input directory
input_dir = inductiva.utils.files.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "reef3d-input-example.zip", unzip=True)
# Initialize the Simulator
reef3d = inductiva.simulators.REEF3D()
# Run simulation with config files in the input directory
task = reef3d.run(input_dir=input_dir, on=mg)
task.wait()
# Once your simulations are done, terminate the machines
mg.terminate()

```
### List active machine groups

You can also list your active machine groups, or get a list of `MachineGroup` objects of previously created machine groups:

```python
import inductiva

inductiva.resources.machine_groups.list()
#                                     Name         VM Type   # machines    Disk Size in GB       Spot         Started at
# api-1b1f724c-5cfe-4d87-8439-9689aa139723   c2-standard-4            1                 40      False   13 Sep, 07:38:50
# api-8e6bf7d8-4888-4de9-bda5-268484b46e6f   c2-standard-4            1                 40      False   13 Sep, 07:37:49

# Get a list of all the MachineGroup objects (for example, if you want to terminate them all at once)
mg_list = inductiva.resources.machine_groups.get()
mg_list
#[<inductiva.resources.machines.MachineGroup at 0x7f8cde53d2a0>,
# <inductiva.resources.machines.MachineGroup at 0x7f8c58954c70>
```
