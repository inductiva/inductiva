# Computational Infrastructure

The Inductiva API serves as a direct intermediary, bridging the gap between users 
and the complex landscape of computational resources. It streamlines the process 
of managing and allocating computational workloads and simulation tasks, by 
facilitating access to an expansive selection of computing options.

This guide will detail how the API simplifies the orchestration of your simulations 
and introduce you to the various computational options currently available to you!

## Computational Resource Management

The Inductiva API acts as an abstraction layer that enables you to access a wide 
array of computational resources provided by a number of different players through 
a unified Python code. These resources could be from cloud providers, bare-metal 
hardware rentals, standard high-performance computing (HPC) solutions commonly 
used in academia, or even on-premise hardware for those with their own computing 
infrastructure. The API serves as a unifying interface atop all these varied resources, 
facilitating access to computational solutions of varying scales, prices, and 
performance levels and helping you select the optimal resource for your needs, 
all through straightforward Python scripting from your laptop.

From the server side, Inductiva manages your computational workload — be it one 
or several simulations. It allocates this workload to the appropriate computational 
resource dedicated to you, handles the orchestration of the simulation, and then 
returns the results back to you.

## Available Computational Resources

Inductiva supports dispatching computational 
workloads to the Google Cloud Platform (GCP). This means that the simulations initiated 
through our API are executed on one or more virtual machines (VMs) hosted on GCP.

There are several families of Virtual Machines (VMs) made available by Inductiva on 
Google Cloud Platform (GCP):

- [**Compute-optimized Machines**]()
- [**Memory-optimized Machines**]()
- [**General-purpose Machines**]()
- [**Accelerator-optimized Machines**]()

To know with more detail which specific machine types Inductiva offers wihtin each machine family, go here.

The VMs within these families are categorized into three types, based on the RAM-to-vCPU 
ratio:

- **highcpu** -  Offers 2 GB of RAM per vCPU, suited for CPU-intensive tasks.
- **standard** -  Provides a balanced 4 GB of RAM per vCPU, ideal for general-purpose use.
- **highmem** - Equipped with 8 GB of RAM per vCPU, designed for memory-intensive applications.

These configurations allow for the customization of
[MachineGroups](http://docs.inductiva.ai/en/latest/api_reference/computational_resources/machinegroup_class.html),
[ElasticMachineGroups](http://docs.inductiva.ai/en/latest/api_reference/computational_resources/elasticgroup_class.html),
and [MPIClusters](http://docs.inductiva.ai/en/latest/api_reference/computational_resources/mpicluster_class.html)
to match your computational needs.

Here's an example of how you can start a MachineGroup with robust "**c3d-standard-60**" 
general machines:

```python
import inductiva

# Initialize a MachineGroup with two "c3d-standard-60" machines
machine = inductiva.resources.MachineGroup(
    machine_type="c3d-standard-60",
    num_machines=2,
)

# Start the MachineGroup
machine.start()
```
Naturally, the cost associated with each machine type varies, and it’s possible 
to accrue significant expenses if a large number of VMs are initiated! To protect 
you from inadvertently spinning up too many resources, the API imposes certain 
limitations on the quantity and types of machines that you can launch. For details 
on these limitations, please consult the
[User Quotas](../basics/quotas.md).

````{eval-rst}
.. seealso::
   Learn how to manage your computational resources through
   `Inductiva's Command Line Interface <../../api-functions/cli/resources.md>`_
````  

## What Next? 

[In our blog post](https://inductiva.ai/blog/article/allocating-computational-resources-in-a-diverse-chip-ecosystem), 
we talked about how the growing diversity of computing options 
is transforming the landscape—and how it’s not always easy to find the best machine 
for your job.

With our [benchmarking tool](benchmarking.md), we’re making it easier
to make smarter, cost-effective decisions for your workloads.