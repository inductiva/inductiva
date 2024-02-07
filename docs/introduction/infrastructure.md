# Computational Infrastructure

## Introduction
One way of understanding the Inductiva API is to see it as an abstraction layer that allows you to use computational resources made available by different providers in a seamless and uniform fashion. 

On the server side, Inductiva takes a computational load issued by the user, i.e. one or several simulations, and routes it to a computational resource that the users has access to, and takes care of all the orchestration required for running the load and returning the results back to the user. 

The computational resources effectively used to process that load can potentially be provided by a number of different players, and it is Inductiva’s role to facilitate the access to those resources, and even help you find the best one for your particular case (e.g. by price or by performance). 

Inductiva API is being designed to allow using compute resources such as:

fully working VMs provided by standard Cloud Providers;
hardware rented by “bare metal” providers;
standard HPC solutions (typically available to academic users);
on-premise hardware (for those who invested in owning their own compute solutions).

Inductiva will function as a uniformizing layer on top of these heterogeneous resources from different providers, allowing users to access computational solutions with different levels of scale, price, performance and flexibility, just by running simple python scripts from their laptop.

## Available Compute Options
As of version 0.4, Inductiva is able to send computational loads to Google Cloud Platform (GCP). Therefore, the simulations you start using the API will run on one VM or several VMs (if you start an MPICluster) that live on GCP. We make available VMs of the following families:

**Compute Optimized:**

c2 https://cloud.google.com/compute/docs/compute-optimized-machines#c2_machine_types
c2d https://cloud.google.com/compute/docs/compute-optimized-machines#c2d_series
h3 https://cloud.google.com/compute/docs/compute-optimized-machines#h3_series

**General Purpose:**

 - [N1]( https://cloud.google.com/compute/docs/general-purpose-machines#n1_machines)
 - N2(https://cloud.google.com/compute/docs/general-purpose-machines#n2_series)
 - N2d - https://cloud.google.com/compute/docs/general-purpose-machines#n2d_machines
 - C3 - 
https://cloud.google.com/compute/docs/general-purpose-machines#c3_series  
- C3d - https://cloud.google.com/compute/docs/general-purpose-machines#c3d_series

Virtual Machines from these families typically come in three variants, with increasing amounts of RAM per vCPU:

- highcpu - 2 Gb per vCPU 
- standard - 4 Gb per vCPU
- highmem - 8 Gb per vCPU

This means you can start MachineGroups, ElasticMachineGroups and MPIClusters with any of these VMs. For example, here is how you start an MachineGroup with beefy c3d general machines:

@ ivan (use standard)

To list the VMs that you can actually request, you can use the CLI:

@ivan 

Of course, each of these machines has a different price, and it is easy to incur in very high-costs if a great number of VM are started. To protect users from accidentally spinning up too many resources, the API imposes certain limitations on how many machines and what type of machines they can launch. Please refer to section [User Quotas]() to know the quotas put in place by the current version of API.

## Other Compute Options
In future versions of the API, we will be making available other compute options that can be selected via the API. Here is a brief description of what is on the roadmap:

Inductiva Compute Engine (ICE). Inductiva is working on its own infrastructure, and will make available cpu-based VMs as well as high-performance GPU-based VMs that certain simulator packages can take advantage of.
Bring your own Cloud. This will allow you to run simulation jobs on your own cloud account (first on GCP account).
Bring your own Hardware. This will allow you to install the API on your own hardware (on-premise or rented) such that you take full advantage of resources you already invested in.
HPC. Inductiva will be working with HPC providers such that certain simulation loads can be redirected to large HPC clusters to which you have access.

These options will be incrementally available to you, over the next releases.

## What to read next

- [User quotas](../user_quotas.md)
