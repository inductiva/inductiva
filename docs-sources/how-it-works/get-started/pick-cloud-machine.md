# Pick the right cloud machine for your simulation

At Inductiva we make available hundreds of cloud machine types.
These machines cover a wide range of specifications, and provide you with options at different levels of performance and cost.

Here’s a short guide to help you navigate this wealth of options.

## Start with the Machine Family

### Leverage your workload

- General-purpose - C4, C4A, C3, C3D, N4, N2, and N2D
Best price-performance ratio and adequate for a variety of workloads.
- Compute-intensive - C2, C2D and H3
Highest performance per core and optimized for compute-intensive workloads.
- Accelerator-optimized - G2, A2 and A3
For workloads that require GPUs.

### Narrow down

- **Generation** - The number identifies the year in which that family was introduced, namely the number “4” shows that the families C4, C4A, N4 were introduced in 2024. Families C3, C3D, H3, A3 are from 2023; and N2, N2D from 2020.
- **Architecture**
  - AMD CPUs
Identified by the letter “D” at the end - e.g. C2D, C3D, N2D.
  - ARM CPUs
End with the letter “A” - e.g. C4A.
  - Equipped with GPUs
Start with “G” for GPU or “A” for Accelerator - e.g. G2, A3, A4.
  - Intel processors - All others, including the ones with GPUs - e.g. C2, N2, C3, H3, C4, G2, A3 and A4.

Once you commit to a certain architecture and key internal specification, it becomes clear what’s the best Machine Family for your case.

---
There are still a few parameters you need to define a specific Machine Type that you can actually use to run your simulations:

- **Number of virtual CPUs (vCPUs)** - More vCPUs allow you to parallelize your simulation more.
- **Number of GPUs** that the machine supports, if it’s the case.
- **Memory profile** or variant of the machine (memory per vCPU) - More RAM per vCPU allows it to handle memory-intensive workloads.
  - highcpu - 1 to 3 GB of RAM per vCPU (typically 2GB)
  - standard - 3 to 7 GB of RAM per vCPU (typically 4GB)
  - highmem - to 14 GB of RAM per vCPU (typically 8GB)
  - megamem - 14 to 19 GB of RAM per vCPU
  - ultramem - 24 to 31 GB of RAM per vCPU

## Leverage RAM and vCPUs

To understand how Inductiva can help you scale your simulations, let’s compare the specs of the machines we make available with a familiar case, such as a good laptop with 32GB of RAM and a 16 core CPU, similar to one that you may own.

At a first glance, such a laptop seems like a reasonably powerful machine and you might be using a similar machine to run relatively large simulations in your office, although sometimes you may have to leave it crunching number for a couple of days.

📈To go from a couple of days to just a few hours, consider the following option that is readily available for you at Inductiva: a c3d-highmem-360 machine. This machine comes equipped with 360 vCPUs (vs the 16 of the laptop) and 8GB of RAM per vCPU, totalling a whopping 2880 GB of RAM. That is 90 more RAM than what is available in the laptop! Very likely, this machine will allow you to run your simulations 5 to 10x faster than in your laptop, and it will, for sure, let you tackle much larger simulation use cases due to its RAM size.
This machine is available for you for about 24,27 US\$ per hour - a good deal for such readily-available compute power - but you can do the job more cost-effectively.

💰If you still need all the CPU firepower but you do not need so much RAM, you can opt for the “highcpu” memory profile, which has less memory per vCPU than the “highmem” profile and, therefore, allows for significant savings. The c3d-highcpu-360 has “only” 720GB of RAM – still more than 20 times the power of the good laptop – but it is 40% less expensive than the previous c3d-highmem-360 option.

---
💡**Inductiva’s Pro tips:**

💰 You can take advantage of the **Spot mode** - Spot VMs have significant discounts, making them great cost savings options for most cases. There is a potential drawback: the cloud provider might stop (preempt) a spot machine at any time to reclaim the compute capacity for other users who are willing to pay the full price. But, if you are running relatively short simulations (e.g. up to 12 hours) then the probability of your spot machine suffering a preemption during that time is very low, making spot instances a great option for saving compute costs on simulation jobs that you can afford to relaunch if needed.

📈 If you privilege raw **performance**, pick a latest-generation machine such as a c4-standard-96, one of the fastest machines we make available, for a regular price of ~5.22 US\$/hour and spot price of ~1.77 US\$/hour).
Newer Generations (e.g. c4, c3, c3d) use more recent technology so they tend to be faster, but they are also more expensive.

💰📈**Top pick**: When in doubt, use the c2d-highcpu-112 machine: it has 112 vCPUs and 224GB of RAM, providing an impressive cost/performance ratio (price around 4.62 US\$/hour and spot price of just ~0.68 US\$/hour).

After narrowing down the machine that best fits your needs, let’s move on and get it started! - <a href="start-first-machine.html">Start your first cloud machine with Inductiva</a>
