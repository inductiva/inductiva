# Simulation Processing and Computational Resource Allocation
In this guide, we will explain some of the main features of the Inductiva API when it comes to managing and utilizing computational resources.

### What We'll Cover
* [Overview]()
* [How can you set up computational resources with Inductiva?]()
    * [Task Queue and Shared VM Pool Integration]()
    * [Custom Hardware Setup for Enhanced Simulation Performance]()
* [MachineGroup Class]()
    * Setting up a MachineGroup for selecting specific VM type
    * Setting up a MachineGroup for running simulations in parallel
* [What to read next]()


## Overview
>*we need something here, I'll think it through - Maya*

## How can you set up computational resources with Inductiva?

>*an intro to the two options - Maya*
### Task Queue and Shared VM Pool Integration

Simulation tasks you create and submit via the Inductiva will be sent to a Task queue, and will eventually be picked up by a pool of workers running on Virtual Machines (VMs) available from a Cloud provider

By default, simulation tasks will be sent to a shared pool of workers serving multiple users. These workers live on VNs that we decided to set aside to make it easier for any user to test the API, and run also to relatively light tasks with the simplest possible setup. 

For example, the code below, will start a SplishSplash simulation that will be automatically picked up by the shared pool of workers:

``````
import inductiva
input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "splishsplash-input-example.zip", unzip=True)


splishsplash = inductiva.simulators.SplishSplash()
task = splishsplash.run(input_dir=input_dir,
                     sim_config_filename="config.json")
task.wait()
task.download_outputs()

``````

>TODO @ivan

>*Observe that at no point we had to explicitly define the target VMs where this simulation would be executed, or even just their specs. Instead, the task will get automatically sent to the shared pool of workers that we prepared for all users. This is very simple, and a great way for doing quick experimentation.*

However, despite being very convenient and simple to use, this pool of resources has a limited predefined capacity, and because it is shared by all users, it is not appropriate for executing larger tasks, since waiting times can be extremely large. So, if you need to run a larger number of simulation tasks, and you need more powerful VMs to run it, you will need to reserve that capacity for your exclusive use.

### Custom Hardware Setup for Enhanced Simulation Performance

Inductiva provides a way of creating pools of VMs resources that are exclusively available for you, and not shared with any other user. We call these Machine Groups, and, as we will explain later, they come in two flavors. 

A Machine Group are groups of homogeneous cloud VMs with specific specs that you can define programmatically and terminate on demand via the API. This will give you full control of the type of VM you use to run your simulations, and will ensure a certain amount of compute power that we reserve exclusively for you. 

Note that a Machine Groups is literally a group of individual VMs that do not communicate with each other. In other words, a MachineGroup is not a cluster, such as an MPI Cluster, where the load of each simulation is divided over all machines of the cluster. To set up an MPI Cluster, see here (point to the MPI documents). 

But for now, let’s dig deeper in the MachineGroup class that implements this notion of a group of independent machines of the same type. 
 
## `MachineGroup` Class

>*@Ivan please explain here the class, each field of the constructor and the main methods etc. Below we will be providing some examples, so here we can be pretty focused on the MachineGroup class.*

Now that we explained what this class is about, let’s see how to use it for two different purposes. First, we will create a MachineGroup with 1 machine only. The goal is to use this functionality to select specific types of VM available on Google Cloud for running our simulation. We will show the impact on performance (and potential cost) of running the simulation on different types of VMs.

Then, we will create a MachineGroup with 5 instances and show how to run 5 variations of the same simulation in parallel. 

### Setting up a Machine Group for selecting specific VM type
>@ivan can we expand the existing example and run it on two different machines.

https://github.com/inductiva/inductiva/blob/development/docs/Machines.md#examples

#### Example

### Setting up a MachineGroup for running simulations in parallel
>@ivan: show a “for loop.”

#### Example

## What to read next
* [Set up an MPI]()
* [Get an overview of the CLI]()
