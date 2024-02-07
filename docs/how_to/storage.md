# Manage Storage

The Inductiva API provides you with tools to manage your remote storage effectively. 
With the Inductiva storage module, you can easily navigate your storage, evaluate the space used, and delete specific directories as needed.
You have the ability to organize your storage by specifying the directory where your simulation outputs should be saved. The directory containing these outputs is automatically named after the task ID. Check [tasks](https://github.com/inductiva/inductiva/tree/main/inductiva/tasks) for more information about this.
Let's illustrate this with some examples.

## Determining the amount of storage in use

As you start running simulations, your remote personal storage will get filled up.
Hence, at times it will be useful to monitor the amount of storage space being utilized.
This can be achieved as follows:

#### Python
```python
import inductiva
space_used = inductiva.storage.get_space_used()
```
#### CLI
```bash
$ inductiva storage size
```

In both cases, we receive the following message:

```
Total user's remote storage in use: 32.4 GB
```

## Viewing storage contents

After determining the total storage space used, you may want to identify which directories are consuming the most storage. To be more specific, you may list
the contents of your directory and sort them by size or creation date as follows:

**Python**

```python
import inductiva
inductiva.storage.listdir(max_results=10, order_by="size", sort_order="desc")
```

**CLI**

```bash
inductiva storage ls --max-results 10 --order-by size --sort-order desc
```

In this way, we obtain a listing of the 10 largest directories within our user's
remote storage.
```bash

       NAME                       SIZE            CREATION TIME
       1699461562982775346/       711.82 MB       08 Nov, 16:39:23
       1699461561456065056/       706.74 MB       08 Nov, 16:39:21
       1698772828190838710/       706.44 MB       31 Oct, 17:20:28
       1698772827312137069/       706.18 MB       31 Oct, 17:20:27
       1699461562181988760/       699.30 MB       08 Nov, 16:39:22
       1698772829464378916/       698.80 MB       31 Oct, 17:20:29
       1699461560739645733/       693.81 MB       08 Nov, 16:39:20
       1698772828793721336/       688.82 MB       31 Oct, 17:20:28
       1698751109500351563/       685.33 MB       31 Oct, 11:18:29
       1699461560029476871/       673.81 MB       08 Nov, 16:39:20
```

In case, you want to be more specific and examine the contents of a specific folder,
you can pass a path to the `listdir` method and/or the CLI subcommand as follows:

**Python**

```python
import inductiva
inductiva.storage.listdir(path = "1234", max_results=10, order_by="size", sort_order="desc")
```

**CLI**

```bash
inductiva storage ls 1699461562982775346/ --max-results 10 --order-by size --sort-order desc
```

Add HERE COMMAND

## Saving simulation outputs

When running a simulation, users can specify the name of the directory where the simulation outputs will be saved on the user's remote storage. This is as simple
as setting the argument `storage_dir` in the `run` method of the simulator. 

This allows users to organize their data on their bucket as they wish.

Let's see an example with the REEF3D simulator.

```python
import inductiva

input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "reef3d-dambreak-example.zip", unzip=True)

simulator = inductiva.simulators.REEF3D()

# Launch the simulation and point the results to the directory "reef3d_simulation"
task = simulator.run(input_dir=input_dir, storage_dir="reef3d_simulation")

task.wait()
task.storage.listdir(order_by="creation_time")
```

```bash
RESULTS
```


### Removing directories

Whenever space needs to be freed up, or in general, when you want to remove a directory
from your remote storage, you can use a simple interface to remove it. In particular,
in case you want to remove everything that is also possible.

The table above provides valuable information that can guide your decision to remove certain directories. 

**Python**

```python
import inductiva
inductiva.storage.rmdir(path="1699461562982775346")
```

**CLI**

```bash
inductiva storage rm 1699461562982775346
```

To remove everything, use the flag `--all`.

The above examples remove the directory `1699461562982775346` permantenly from the user's remote storage. So be careful when using this command.


### Machine Groups
You can control the amount of VM storage dedicated to storing the results of your 
simulation using the parameter `data_disk_gb` of the MachineGroup class. You set
this parameter when you instantiate the MachineGroup object, and this becomes fixed
for the corresponding VMs since it is not possible to change storage allocation
after instantiation. Below is an example of how you would reserve 20GB of storage
in each machine when starting a MachineGroup with 5 machines:

```python
import inductiva

machine = inductiva.resources.MachineGroup(
    "c2-standard-16", num_machines=5, data_disk_gb=20)
```

### MPI Clusters

For MPI Clusters, the machines in the cluster share an NFS partition where simulators 
typically write their final results. So, in this case, the storage parameter `data_disk_gb`
sets the storage size for the NFS partition. Here is an example where we set an MPI 
Cluster with 8 machines, where each machine has been given 50Gb of NFS storage:

```python
import inductiva

mpi_cluster = inductiva.resources.MPICluster(
    "c2-standard-16", num_machines=8, data_disk_gb=50)
```

### What about machines in the Common Pool?

Machines in the Common Pool have a storage space of 30 GB that you can’t control. This 
means that if you are submitting to the Common Pool simulations that produce more than 
30 GB, they will fail, and your Task will fail. There is no way for you to request more 
storage space for VMs in the Common Pool. Common Pool machines are intended for running 
short simulations, mostly with the goal of testing your scripts. If you wish to run 
simulations that produce a large amount of data, then you really need to spin up your own 
Machine Groups or MPI Clusters.


