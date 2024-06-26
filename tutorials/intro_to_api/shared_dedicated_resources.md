# Shared and Dedicated Resources

In this guide, we will explain some of the main features of the Inductiva API when 
it comes to making informed decisions about resource allocation. Here, you will 
learn the two options available for running your simulations; shared and dedicated
resources. Additionally, we will demonstrate the performance differences between 
the two resources using a SWASH simulation.

## Shared Resources

By default, when you submit simulation tasks via the Inductiva API, they are queued 
and dispatched to a shared pool of workers on designated Virtual Machines (VMs) 
serving multiple users, primarily utilizing resources from the
[Google Cloud Provider (GCP)](https://cloud.google.com/compute/docs/machine-resource). 
Looking ahead, future versions will also support Inductiva's own computational platform.

These VMs are specifically allocated to facilitate easy API testing and to handle light 
tasks efficiently with a minimal setup, ideal for quick experimentation.

However, this shared resource pool has its limitations, including slower completion 
times for simulations due to its finite capacity and less powerful VMs. For running 
a larger volume of simulations, the Inductiva API provides you with the option to 
set up your own dedicated resources. 

## Dedicated Resources

To ensure you can access the necessary computational power for your simulations 
without the constraints of a shared environment, Inductiva provides the option 
to create dedicated pools of Virtual Machine (VM) resources, exclusively 
reserved for your use and not shared with others. 

In this option, there are three types of dedicated computational resources you can
launch for your simulations:

- [**Machine Group**](https://docs.inductiva.ai/en/latest/api_reference/computational_resources/machinegroup_class.html): 
This consists of homogeneous machines designed to operate individually, enabling 
the distribution of multiple simulations across different machines for parallel processing.
- [**Elastic Machine Group**](http://docs.inductiva.ai/en/latest/api_reference/computational_resources/elasticgroup_class.html): 
Similar to Machine Group, these also consist of individual machines. The key difference here is the 
elastic scaling feature, which dynamically adjusts the number of machines based 
on simulation demands, ensuring efficient resource utilization.
- [**MPI Cluster**](http://docs.inductiva.ai/en/latest/api_reference/computational_resources/mpicluster_class.html) 
This setup involves a network of machines configured to work in tandem on a single 
simulation task, distributing the workload across multiple CPUs. This is particularly 
useful for complex simulations that exceed the capabilities of a single machine.

## Comparison Example: SWASH

To illustrate the performance differences between running simulations on a shared 
pool of resources and on dedicated resources using the Inductiva API, we will run 
[SWASH](https://docs.inductiva.ai/en/latest/simulators/SWASH.html) as an example.

### SWASH on Shared Resources

First, we'll run the simulation using the shared pool of workers, a convenient 
option for those getting started or running less resource-intensive tasks. 

```python
import inductiva

# Download the input files for the SWASH simulation
input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "swash-resources-example.zip", unzip=True)

# Initialize the SWASH simulator and run the simulation
# in the default shared pool of workers
swash = inductiva.simulators.SWASH()
task = swash.run(input_dir=input_dir,
                 sim_config_filename="input.sws")

# Wait until the task finishes running.
task.wait()
```
In this option, it takes **approximately 20 minutes** to complete this simulation.

### SWASH on Dedicated Resources

Let's now run the same simulation on dedicated resources, specifically set 
up for this task.

```python
import inductiva

# Download the input files for the SWASH simulation
input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "swash-resources-example.zip", unzip=True)

# Instantiate a MachineGroup object with 1 preemptible machine of type
# c2-standard-30 and start it immediately
machine_group = inductiva.resources.MachineGroup(
    machine_type="c2-standard-30", spot=True)
machine_group.start()

# Initialize the SWASH simulator and run the simulation
# in your just launched dedicated MachineGroup
swash = inductiva.simulators.SWASH()
task = swash.run(input_dir=input_dir,
                 sim_config_filename="input.sws",
                 on=machine_group)

# Wait for the task to finish and download the outputs
task.wait()

# Terminate your dedicated MachineGroup at then end of the simulation.
machine_group.terminate()
```
Notice that, the simulation is picked almost immediately - no waiting time required,
and by contrasting the two resource options above, running the simulation on a
[dedicated machine group](#dedicated-resources) with a `c2-standard-30` machine
took **9 minutes and 37s** to complete, which is **2.68 times less** than on the 
[shared pool](#shared-resources). 

