# Simulation Processing and Computational Resource Allocation
In this guide, we will explain some of the main features of the Inductiva API when it comes to managing and utilizing computational resources.You will learn how to set up computational resources that will be exclusively dedicated to running your simulation tasks, and that will allow you to scale up your projects to hundreds of machines.

### What We'll Cover
* [Overview]()
* [How can you set up computational resources with Inductiva?]()
    * [Task Queue and Shared VM Pool Integration]()
    * [Custom Hardware Setup for Enhanced Simulation Performance]()
* [MachineGroup Class]()
    * Setting up a MachineGroup for selecting specific VM type
    * Setting up a MachineGroup for running simulations in parallel
* [What to Read Next]()


## Overview


## How can you set up computational resources with Inductiva?

### Task Queue and Shared VM Pool Integration
By default, simulation requests are processed by a shared pool of machines serving 
multiple users. These machines are equipped with 4 cores each, and the pool is 
limited to a predefined capacity. 

Simulation tasks you create and submit via the Inductiva will be sent to a Task queue, and will eventually be picked up by a pool of workers running on Virtual Machines (VMs) available from a Cloud provider. 

By default, simulation tasks will be sent to a shared pool of workers serving multiple users. These workers live on VNs that we decided to set aside to make it easier for any user to test the API, and run also to relatively light tasks with the simplest possible setup. 

### Custom Hardware Setup for Enhanced Simulation Performance
If you require dedicated resources with custom 
configurations, the **Inductiva API** allows you to easily set up virtual machines 
reserved solely for your simulations. These are managed via machine groups, *i.e.*, 
groups of homogeneous machines with specific requested properties that can be 
started and terminated on demand. 

Creating a `MachineGroup` requires specifying the 
type of Virtual Machine to use and the number of such VMs that will compose the 
group. Currently, the available options for the `machine_type` are the ones 
available in the [Google Cloud Platform](https://cloud.google.com/compute/docs/machine-types). 

Once a `MachineGroup` is created, simply pass it as argument to your simulations, 
which will then be scheduled to run on those machines. 

## MachineGroup Class

With the `MachineGroup` class, you can configure, 
start, and terminate machines. 

choose one of the options below

### Setting up a Machine Group for selecting specific VM type

#### Example

### Setting up a MachineGroup for running simulations in parallel

#### Example

