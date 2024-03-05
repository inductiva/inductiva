# Computational Infrastructure

## Introduction
One way of understanding the Inductiva API is to see it as an abstraction layer that
allows you to use computational resources made available by different providers in a
seamless and uniform fashion. 

On the server side, Inductiva takes a computational load issued by the user, i.e.,
one or several simulations, and routes it to a computational resource that the users
has access to, and takes care of all the orchestration required for running the load
and returning the results back to the user. 

The computational resources effectively used to process that load can potentially be
provided by a number of different players, and it is Inductiva’s role to facilitate
the access to those resources, and even help you find the best one for your particular
case (e.g. by price or by performance). 

Inductiva API is being designed to allow using compute resources such as:

- fully working VMs provided by standard Cloud Providers;
- hardware rented by “bare metal” providers;
- standard HPC solutions (typically available to academic users);
- on-premise hardware (for those who invested in owning their own compute solutions).

Inductiva will function as a uniformizing layer on top of these heterogeneous resources
from different providers, allowing users to access computational solutions with different
levels of scale, price, performance and flexibility, just by running simple python
scripts from their laptop.

## Available Computational Options
As of version 0.4, Inductiva is able to send computational loads to Google Cloud
Platform (GCP). Therefore, the simulations you start using the API will run on one
VM or several VMs (if you start an MPICluster) that live on GCP. We make available
VMs of the following families:

**Compute Optimized:**

- [C2](https://cloud.google.com/compute/docs/compute-optimized-machines#c2_machine_types)
- [C2D](https://cloud.google.com/compute/docs/compute-optimized-machines#c2d_series)
- [H3](https://cloud.google.com/compute/docs/compute-optimized-machines#h3_series)

**General Purpose:**

 - [N1]( https://cloud.google.com/compute/docs/general-purpose-machines#n1_machines)
 - [N2](https://cloud.google.com/compute/docs/general-purpose-machines#n2_series)
 - [N2D](https://cloud.google.com/compute/docs/general-purpose-machines#n2d_machines)
 - [C3](https://cloud.google.com/compute/docs/general-purpose-machines#c3_series)  
 - [C3D](https://cloud.google.com/compute/docs/general-purpose-machines#c3d_series)

Virtual Machines from these families typically come in three variants, with increasing amounts of RAM per vCPU:

- **highcpu** - 2 Gb per vCPU 
- **standard** - 4 Gb per vCPU
- **highmem** - 8 Gb per vCPU

This means you can start MachineGroups, ElasticMachineGroups and MPIClusters with any of these VMs. For example, here is how you start a MachineGroup with beefy C3D general machines:

```python
import inductiva

machine = inductiva.resources.MachineGroup(
    machine_type="c3d-standard-60",
    num_machines=2,
)

machine.start()
```

To list the VMs that you can actually request, you can use the CLI:

```
$ inductiva resources available
Machine types provided in Google Cloud

c2: Intel Xeon Cascade Lake (2nd Gen) processor.
  > c2-standard-  [2, 4, 8, 16, 30, 60]                         

c3: Intel Xeon Sapphire Rapids (4th Gen) processor.
  > c3-highcpu-   [4, 8, 22, 44, 88, 176]                       
  > c3-standard-  [4, 8, 22, 44, 88, 176]                       
  > c3-highmem-   [4, 8, 22, 44, 88, 176]                       

h3: (Available Soon) Intel Xeon Sapphire Rapids (4th Gen) processor.
Simultaneous multithreading disabled, i.e., vCPU represents an entire core.
  > h3-standard-  [88]                                          

c2d: AMD EPYC Milan (3rd Gen) processor.
  > c2d-highcpu-  [2, 4, 8, 16, 32, 56, 112]                    
  > c2d-standard- [2, 4, 8, 16, 32, 56, 112]                    
  > c2d-highmem-  [2, 4, 8, 16, 32, 56, 112]                    

c3d: AMD EPYC Genoa (4th Gen) processor.
  > c3d-highcpu-  [4, 8, 16, 30, 60, 90, 180, 360]              
  > c3d-standard- [4, 8, 16, 30, 60, 90, 180, 360]              
  > c3d-highmem-  [4, 8, 16, 30, 60, 90, 180, 360]              

e2: Intel Xeon (up to Skylake, 1st Gen) and AMD EPYC (up to Milan, 3rd Gen) processors.
Automatically selected based on availability.
  > e2-highcpu-   [2, 4, 8, 16, 32]                             
  > e2-standard-  [2, 4, 8, 16, 32]                             
  > e2-highmem-   [2, 4, 8, 16]                                 

n2: Intel Xeon Ice Lake and Cascade Lake processors (3rd and 2nd Gen).
Cascade Lake default up to 80 vCPUs and Ice Lake for larger machines.
  > n2-highcpu-   [2, 4, 8, 16, 32, 48, 64, 80, 96]             
  > n2-standard-  [2, 4, 8, 16, 32, 48, 64, 80, 96, 128]        
  > n2-highmem-   [2, 4, 8, 16, 32, 48, 64, 80, 96, 128]        

n2d: AMD EPYC Milan or ROME processors (3rd and 2nd Gen).
  > n2d-highcpu-  [2, 4, 8, 16, 32, 48, 64, 80, 96, 128, 224]   
  > n2d-standard- [2, 4, 8, 16, 32, 48, 64, 80, 96, 128, 224]   
  > n2d-highmem-  [2, 4, 8, 16, 32, 48, 64, 80, 96]             

n1: Intel Xeon (up to Skylake, 1st Gen) processor.
Automatically selected based on availability.
  > n1-highcpu-   [1, 2, 4, 8, 16, 32, 64, 96]                  
  > n1-standard-  [1, 2, 4, 8, 16, 32, 64, 96]                  
  > n1-highmem-   [1, 2, 4, 8, 16, 32, 64, 96]                  
```

In case, you more information you can use the `-v` flag, and to focus on a specific series, you can use:

```
$ inductiva resources available -s c3d -v
Machine types provided in Google Cloud

c3d: AMD EPYC Genoa (4th Gen) processor.
  > c3d-highcpu-  [4, 8, 16, 30, 60, 90, 180, 360]              -> 2 GB of memory per vCPU.
  > c3d-standard- [4, 8, 16, 30, 60, 90, 180, 360]              -> 4 GB of memory per vCPU and possible local ssd integration.
  > c3d-highmem-  [4, 8, 16, 30, 60, 90, 180, 360]              -> 4 GB of memory per vCPU.

```

Of course, each of these machines has a different price, and it is easy to 
incur in very high costs if a great number of VM are started. To protect users 
from accidentally spinning up too many resources, the API imposes certain 
limitations on how many machines and what type of machines they can launch. 
Please refer to the section [User Quotas](../api_reference/user_quotas.md) to 
know the quotas put in place by the current version of API.

## Other Computational Options
In future versions of the API, we will be making available other compute options
that can be selected via the API. Here is a brief description of what is on the
roadmap:

- **Inductiva Compute Engine (ICE)**. Inductiva is working on its infrastructure
and will make available cpu-based VMs as well as high-performance GPU-based VMs
that certain simulator packages can take advantage of.
- **Bring your own Cloud (BYOC)**. This will allow you to run simulation jobs on
your own cloud account (first on GCP account).
- **Bring your own Hardware (BYOH)**. This will allow you to install the API on
your own hardware (on-premise or rented) such that you take full advantage of
resources you already invested in.
- **High-performance Computing (HPC)**. Inductiva will be working with HPC providers
such that certain simulation loads can be redirected to large HPC clusters to which
you have access.

These options will be incrementally available to you, over the next releases.

## What to read next

- [User Quotas](../api_reference/user_quotas.md)
