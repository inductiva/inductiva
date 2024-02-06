
> TO BE BROKEN INTO PIECES

# Computational Resources

Besides running and managing simulations, Inductiva API enhances the computational 
infrastructure of any user by providing a simple interface to run simulations in
your own dedicated computational resources.

In this tutorial, we learn the API capabilities to launch and manage 
computational resources in the cloud and guide the scaling of your 
computational infrastructure to run simulations faster and more efficiently.

Before getting into the details, we highlight that launching a computational 
resource and running a simulation in it is as simple as follows:

```python
import inductiva

# Configure a computational resources
machines = inductiva.resources.MachineGroup(
    machine_type="c2-standard-16", num_machines=5, data_disk_gb=60)

# Launch the computational resource to be available to run simulations
machines.start()

# Run the Getting Started example simulation in our own dedicated computational
# resource
input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "reef3d-dambreak-dir.zip", unzip=True)

simulator = inductiva.simulators.REEF3D()

# Select the computational resource to run the simulation on the `run` method
task = simulator.run(input_dir=input_dir, on=machines)

task.wait()

# Terminate the computational resource when you don't need it anymore
machines.terminate()
```

This code snippet launches a computational resource with 5 machines of type
`c2-standard-16` that are available to now deploy your own simulations. In the
following sections, we will dive deep into the details of the computational resources
available via the API and how to use them.

### Steps we will cover:
1. [Set up a MachineGroup for running simulations in parallel]()
2. [Launch a Machine Group]()
2. [Set up an Elastic Machine Group]()
4. [Start a MPI Cluster in the Cloud]()
5. [Manage the computational resources]()


## Set up a MachineGroup for running simulations in parallel

#### Example

The second use case of launching a `MachineGroup` is that of setting multiple
simulations running in parallel and distributed by the various machines constituting
the group. This is useful when you want to run multiple simulations in parallel,
but you don't want to wait for the first one to finish before starting the second one. 

To exemplify, we will use the [templating mechanism]() built-in the Inductiva API
to automatically change the water level of the simulation in the input files and
run 5 different simulations in parallel. 

```python
import inductiva
from inductiva import mixins

# Instantiate a MachineGroup object with 1 preemptible machine of type
# c2-standard-30 and start it immediately
machine_group = inductiva.resources.MachineGroup(
    machine_type="c2-standard-30", num_machines=5, spot=True)
machine_group.start()

# Download the input files for the SWASH simulation
input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "swash-template-example.zip", unzip=True)

# Initialize the template file manager
file_manager = mixins.FileManager()

# Initialize the SWASH simulator
swash = inductiva.simulators.SWASH()

# Explore the simulation for different water levels
water_levels_list = [3.5, 3.75, 4.0, 4.5, 5.0]

# Launch multiple simulations
for water_level in water_levels_list:
    # Set the root directory and render the template files into it.
    file_manager.set_root_dir("swash-input-example")
    file_manager.add_dir(input_dir, water_level=water_level)

    # Run the simulation on the dedicated MachineGroup
    task = swash.run(input_dir=file_manager.get_root_dir(),
                    sim_config_filename="input.sws",
                    on=machine_group)
```

The template mechanism will allow you to explore 5 different variations of the
simulation, each with a different water level. The simulations will be submitted
to our dedicated machine group and will run in parallel.
We can check that all simulations are running via the CLI and that it took only
1min for the moment they are submitted until they start running:

```bash
$ inductiva tasks list
ID                         Simulator    Status    Submitted         Started           Computation Time    Resource Type
-------------------------  -----------  --------  ----------------  ----------------  ------------------  ---------------
57mr4kas99jxb9titkeackano  swash        started   01 Feb, 09:07:19  01 Feb, 09:08:03  *0:03:12            c2-standard-30
ox8718m0pwfi02zczui3qky4w  swash        started   01 Feb, 09:07:17  01 Feb, 09:08:02  *0:03:14            c2-standard-30
mak1ji62s7axf7mespkc36g7e  swash        started   01 Feb, 09:07:15  01 Feb, 09:08:03  *0:03:14            c2-standard-30
ijyu8bkvme7vg9k0kj6v23gxa  swash        started   01 Feb, 09:07:14  01 Feb, 09:08:02  *0:03:16            c2-standard-30
g5qq5c9mk2nr5wqhzef38sdm4  swash        started   01 Feb, 09:07:12  01 Feb, 009:08:01  *0:03:17            c2-standard-30
```

This is a great way to speed up the execution of multiple simulations, since the
time to run all 5 simulations will be approximately the same as running just one,
that is the above 5 simulations took 9m55s to complete, which is the time of
the slowest simulation.

Now, that all the simulations have finished running, we end this tutorial with an
extra lesson to help reduce the amount of time machines are left unused:
> Don't forget to terminate your computational resources with `inductiva resources terminate --all`.

#### Configuration parameters

The above computational resources are configured with a few common parameters. Let's
start by introducing them:

- `machine_type` defines the type of CPU used for each machine. This parameter
follows the naming convention set by [Google Cloud](https://cloud.google.com/compute/docs/machine-types).
Each machine type, e.g. `c2-standard-60`, is composed of a prefix defining the
CPU series, a suffix that sets the number of [virtual CPUs (vCPU)](https://cloud.google.com/compute/docs/cpu-platforms)
per machine and the middle word refers to the level of RAM per vCPU. In the example,
`c2` refers to an Intel Xeon Scalable processor of 2nd generation, `standard`
means 4 GB of RAM per vCPU and will contain `60` vCPUs. Below, we introduce the
machine types available via the API. 
- the parameters `num_machines`, `min_machines`, `max_machines` control the number
of machines available in the computational resource. The former is used by the
Machine group and the MPI Cluster resources. The latter ones set the minimum number
of machines that are always available and the maximum number of machines that can
be available at any time, respectively. 
- `data_disk_gb` allows the selection of the size of the disk attached to each machine that is reserved for the simulation data in GB.
- `spot` allows to differ between standard resources or preemptible ones. The latter
resources are way cheaper but only serve for fault-tolerant workloads since they
can be stopped at any time. The former are fully dedicated to the user's usage.
In case of failure of a simulation with a preemptible machine, the simulation is
resubmitted to the queue of the computational resource. Currently, this is only
available for the Machine Group and Elastic Machine Group.


## Launch a Machine Group

The first computational resource we introduce is the `MachineGroup`, which are a
group of homogeneous machines that do not communicate with each other and work
individually to run multiple simulations simultaneously. The number of machines
in the group is fixed and the machines are always active to run simulations.
Based on the parameters above you can configure it as you wish, and then start
it to make it available for your simulations.

Launching a `MachineGroup` is as simple as follows:

```python
import inductiva

# Configure a simple machine group with 2 preemptibles machines, each with a disk
# size of 60 GB and an Intel Xeon Scalable processor of 2nd gen with 8 vCPUs and
# 4 GB of RAM per vCPU.
machine_group = inductiva.resources.MachineGroup(
    machine_type="c2-standard-8",
    num_machines=2,
    data_disk_gb=60,
    spot=True
)

# Launch the machine group to make it available to run simulations:
machine_group.start()
```

Once a `MachineGroup` is created, simply pass it as an argument on the simulator 
`run` method, and it will be added to the queue of the computational resource.


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

## Start an MPI Cluster in the Cloud

To empower your computational infrastructure further, the other computational resource
available is the `MPICluster`. This is a group of homogenous machines that
work together to run a single simulation at a time. An MPI cluster is used when a
simulation requires more computational power than a single machine can provide.
Usually, launching an MPI cluster locally involves several steps to make sure the
machines are communicating with each other to collaborate on the simulation. The
API takes care of all the magic to make sure that everything is set up correctly so
that you can focus on running your simulations.

The `MPICluster` is configured with the same parameters as the `MachineGroup`, but
at the moment the machines can not be preemptable. 

```python
import inductiva

# Configure an MPICluster with 3 machines, each with a disk size of 70 GB and an
# AMD EPYC Milan 3rd gen processor with 16 vCPUs  Intel Xeon Scalable processor of 
# 2nd gen with 16 vCPUs and 8 GB of RAM per vCPU.
mpi_cluster = inductiva.resources.MPICluster(
    machine_type="c2d-highmem-16", num_machines=3, data_disk_gb=70)

# Launch the MPI cluster to make it available to run simulations:
mpi_cluster.start()
```

With the `MPICluster` started you are ready to launch massive simulations to get
your cluster working together.


## Manage the computational resources

Once you have launched your computational resources, there are a few API methods
that help manage them and see their status. Let's go over them one by one.

#### Get your active computational resources

In case computational resources have been launched in another Python session and
want to manage or re-use them in another script, users can fetch the respective
instance via the `get` method as follows:

```python
## Obtain a list with instances of all active computational resources
>>> resources_list = inductiva.resources.machine_groups.get()
>>> print(resources_list)
[MPICluster(name="api-23zssj6oq77xxsot3o0nhax3d"),
 ElasticMachineGroup(name="api-45fetsr58okcs0x6j9m0vsi2z"),
 MachineGroup(name="api-4kken08fnoxuu5zjjak6ak2xe")]
```

#### List your active computational resources

When you just want to check the active resources you can quickly list the
information of each one either via Python or via the CLI.

Via Python
```python
inductiva.resources.machine_groups.list()
```
or via CLI:
```bash
$ inductiva resources list
```

One obtains for example the following information:
```
Active Resources:
Name                           Machine Type    Elastic    Type        # machines    Disk Size in GB  Spot    Started at (UTC)
-----------------------------  --------------  ---------  --------  ------------  -----------------  ------  ------------------
api-23zssj6oq77xxsot3o0nhax3d  c2d-highmem-16  False      mpi                  3                 70  False   01 Feb, 12:30:06
api-45fetsr58okcs0x6j9m0vsi2z  c2-standard-4   True       standard           1/5                 70  False   01 Feb, 12:25:54
api-4kken08fnoxuu5zjjak6ak2xe  c2-standard-8  False      standard              2                 60  True    01 Feb, 12:26:37
```


#### Terminate the active computational resources

When you have finished using your computational resources, don't forget to terminate
them. An advantage of Inductiva API is to control your computational
resources and avoid them being idle.

Hence, you can either terminate your computational resources via Python or the CLI. 

In Python, you will need to instantiate the computational resource object, for
example, with the `get` method and then do `machine.terminate()` for example. 

Via CLI, this process is much easier and you can either terminate a machine on
at a time:
```bash
inductiva resources terminate api-45fetsr58okcs0x6j9m0vsi2z
```

Or terminate them all at once. You will need to confirm this action:
```bash
inductiva resources terminate --all
```

Both are a blocking call that will only finish when the machines have terminated,
in this way no computational resources are left up.
