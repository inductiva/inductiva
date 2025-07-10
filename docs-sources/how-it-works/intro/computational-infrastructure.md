# Computational Infrastructure

The Inductiva API serves as a direct intermediary, bridging the gap between users 
and the complex landscape of computational resources. It streamlines the process 
of managing and allocating computational workloads and simulation tasks, by 
facilitating access to an expansive selection of computing options.

This guide will detail how the API simplifies the orchestration of your simulations 
and introduce you to the various computational options currently available to you!

## Computational Resource Management

The Inductiva API provides a unified Python interface that simplifies access to diverse computational resources across multiple providers. These resources could be from cloud providers, physical server rentals, standard high-performance computing (HPC) solutions commonly used in academia, or even on-premise hardware for those with their own computing infrastructure. The API serves as a unifying interface atop all these varied resources, facilitating access to computational solutions of varying scales, prices, and 
performance levels and helping you select the optimal resource for your needs, 
all through straightforward Python scripting from your laptop.

From the server side, Inductiva manages your computational workload — be it one 
or several simulations. It allocates this workload to the appropriate computational 
resource dedicated to you, handles the orchestration of the simulation, and then 
returns the results back to you.

## Available Computational Resources

By default, Inductiva executes computational workloads on Google Cloud Platform (GCP). This means that the simulations initiated through our API are executed on one or more virtual machines (VMs) hosted on GCP.

There are several families of Virtual Machines (VMs) made available by Inductiva on 
Google Cloud Platform (GCP):

- [**Compute-optimized Machines**](https://cloud.google.com/compute/docs/compute-optimized-machines): Ideal for CPU-intensive simulations requiring high-performance processors and optimized compute-to-memory ratios.
- [**Memory-optimized Machines**](https://cloud.google.com/compute/docs/memory-optimized-machines): Perfect for memory-intensive applications that require large amounts of RAM.
- [**General-purpose Machines**](https://cloud.google.com/compute/docs/general-purpose-machines): Versatile VMs that provide a balanced mix of compute, memory, and networking resources, suitable for a wide range of simulation workloads.
- [**Accelerator-optimized Machines**](https://cloud.google.com/compute/docs/accelerator-optimized-machines): Specialized VMs equipped with GPUs and other accelerators for high-performance computing tasks that benefit from parallel processing capabilities.

Each VM family offers multiple machine types with different specifications to match your simulation requirements and budget. For detailed information about specific machine types available within each family through Inductiva, including pricing and performance characteristics, explore our complete machine catalog.

These machine types allow for the customization of
[MachineGroups](computational_resources/machinegroup_class.md),
[ElasticMachineGroups](computational_resources/elasticgroup_class.md),
and [MPIClusters](computational_resources/mpicluster_class.md)
to match your computational needs.

Here's an example of how you can start a MachineGroup using "c3d-standard-60" machines:

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
   `Inductiva's Command Line Interface <https://inductiva.ai/guides/documentation/cli/resources>`_
```` 

### Beyond Google Could Platform

While the examples above focus on GCP resources, Inductiva also supports running simulations on other computational infrastructures. For comprehensive information about using your own hardware, see our [BYOH (Bring Your Own Hardware)](../../expand/bring-your-own-software/index.md) tutorial, which explores in detail the various infrastructure options we support beyond GCP.

## What Next? 

[In our blog post](https://inductiva.ai/blog/article/allocating-computational-resources-in-a-diverse-chip-ecosystem), 
we talked about how the growing diversity of computing options 
is transforming the landscape—and how it’s not always easy to find the best machine 
for your job.

With our [benchmarking tool](benchmarking.md), we make it easier to make smarter, cost-effective decisions for your workloads.