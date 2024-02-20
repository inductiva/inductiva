# Computational Infrastructure

One way to understand the Inductiva API is to see it as an intermediary between 
users and the complex landscape of computational resources. In this reference, you 
will learn how the Inductiva API manages your computational workload, and then you'll
explore the various computational options currently available to you, and additional resources set to be introduced in future updates.

## Computational Resource Management

The Inductiva API acts as an abstraction layer that enables you to access a wide 
array of computational resources across different providers through a unified 
interface called the [Command Line Interface (CLI)](). The API's key role is to 
assist you in identifying the resources best suited to your needs, whether you're prioritizing cost-efficiency or performance.

These resources could be from cloud providers, bare-metal hardware rentals, standard 
high-performance computing (HPC) solutions commonly used in academia, or even on-premise 
hardware for those with their own computing infrastructure. The API serves as a unifying interface atop all these varied resources, allowing you to tap into computational solutions of varying scales, prices, and performance levels through straightforward Python 
scripting from your laptop.

From the server side, Inductiva manages your computational workload — be it one 
or several simulations. It allocates this workload to the appropriate computational 
resource available to you, handles the orchestration of the simulation, and then 
returns the results back to you.

## Available Computational Resources

In the current release, version 0.4, Inductiva supports dispatching computational 
workloads to the Google Cloud Platform (GCP). This means that the simulations initiated through our API are executed on one or more virtual machines (VMs) hosted on GCP. 

The VM families made available by Inductiva on Google Cloud
Platform (GCP) include:


````{eval-rst}
.. tabs::

   .. tab:: Compute-optimized Machines

      - [C2](https://cloud.google.com/compute/docs/compute-optimized-machines#c2_machine_types)
      - [C2D](https://cloud.google.com/compute/docs/compute-optimized-machines#c2d_series)
      - [H3](https://cloud.google.com/compute/docs/compute-optimized-machines#h3_series)

   .. tab:: General-purpose Machines

       - [N1]( https://cloud.google.com/compute/docs/general-purpose-machines#n1_machines)
       - [N2](https://cloud.google.com/compute/docs/general-purpose-machines#n2_series)
       - [N2D](https://cloud.google.com/compute/docs/general-purpose-machines#n2d_machines)
       - [C3](https://cloud.google.com/compute/docs/general-purpose-machines#c3_series)  
       - [C3D](https://cloud.google.com/compute/docs/general-purpose-machines#c3d_series)

````

The VMs within these families are categorized into three types, based on the RAM-to-vCPU ratio:

- **highcpu** -  Offers 2 GB of RAM per vCPU, suited for CPU-intensive tasks.
- **standard** -  Provides a balanced 4 GB of RAM per vCPU, ideal for general-purpose use.
- **highmem** - Equipped with 8 GB of RAM per vCPU, designed for memory-intensive applications.

These configurations allow for the customization of [MachineGroups](), [ElasticMachineGroups](), and [MPIClusters]() to match your computational needs.

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
Naturally, the cost associated with each machine type varies, and it's possible 
to accrue significant expenses if a large number of VMs are initiated! To protect
you from inadvertently spinning up too many resources, the API imposes certain 
limitations on the quantity and types of machines that you can launch. For details 
on these limitations, please consult the [User Quotas](../api_reference/user_quotas.md) section, to go over the quotas we put in place through the current version of the API.

## Upcoming Computational Resources

In future versions of the API, we plan to expand the range of computational resources 
you can access through the API. Here's a sneak peek at what is on our development roadmap:

- **Inductiva Compute Engine (ICE):** We are developing our own infrastructure to 
offer both CPU-based and high-performance GPU-based VMs. These will be particularly beneficial for simulator packages requiring extensive computational power.

- **Bring your own Cloud (BYOC):** This feature will enable you to use your own cloud accounts, starting with GCP, to run simulation jobs.

- **Bring your own Hardware (BYOH):** With BYOH, you'll be able to integrate the 
API directly with your existing hardware, whether it's on-premise or leased, 
maximizing the resources you already invested in.

- **High-performance Computing (HPC):** We are collaborating with HPC providers 
to allow redirection of certain simulation loads to large HPC clusters you can
access.

Stay tuned for these updates, which will be rolled out in the upcoming releases!


## What to read next

- [User Quotas](../api_reference/user_quotas.md)
