---
myst:
  html_meta:
    description: "Explore available computational resources through the Inductiva API and learn how to choose more powerful hardware to run the 'base case'"
    keywords: "Inductiva API, Programming, HPC, Simulation, Tutorial, Synthetic Data Generation, Physics-ML, SPH"
---

# Choosing your Hardware

## Understanding what is happening under the hood

In simulation projects, efficient access to the right computational resources is
crucial for achieving accurate and timely results. Our API offers users the
flexibility to optimize their workflows by **launching dedicated machines**
tailored to their specific needs. These machines come in a range of configurations,
from high-performance options for intensive simulations to more cost-effective
setups for lighter tasks. This allows users to balance computational power and
budget as required.

By creating a MachineGroup, users can run simulations efficiently, adjusting
parameters like particle count for more complex scenarios, all while ensuring
execution times remain manageable. This adaptability makes it easier to scale
computations as project demands grow.

## Browsing our Available Machines

The best way to speed things up is to use our API to launch more powerful
user-dedicated machines, reserved for our own exclusive use.

Let's explore alternative machine setups available via the API by running a simple command through the Inductiva [Command Line Interface (CLI)](https://docs.inductiva.ai/en/latest/cli/cli-overview.html):

```console
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

e2: Intel Xeon (up to Skylake, 1st Gen) and AMD EPYC (up to Milan, 3rd Gen)
processors.
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

OK! This looks good. Let's choose something different from the menu.

## Re-Running our "Base Case" on a better VM

We'll run our "base case" with a _particle radius of 0.008_ again, but now we'll
boost our computing power by scaling up the number of vCPUs in the `c2` family.
We're moving from the 4 vCPUs of the `c2-standard-4` to a much heftier `c2-standard-60`
setup, which hopefully will allow us to significantly decrease the simulation time.
This can easily be achieved with just a few extra lines of code:

```python
import inductiva

# Configure and start a dedicated machine group
my_machine_group = inductiva.resources.MachineGroup(
    machine_type="c2-standard-60")
my_machine_group.start()

input_dir = "splishsplash-base-dir"

# Set simulation input directory
splishsplash = inductiva.simulators.SplishSplash()

task = splishsplash.run(input_dir=input_dir,
                        sim_config_filename="config.json",
                        on=my_machine_group)

task.wait()

# Terminate the machine group
my_machine_group.terminate()

task.download_outputs()
```

Notice the significant reduction in runtime for this simulation: down from 30m to just **10m07s**, achieving an approximately 3-fold increase in speed! Remember that, in practice, speed ups are not linear with the number of vCPUs, so this is already quite an achievement with just a few lines of code.

We can enhance the performance further by choosing even more powerful machines. For example, the `c3` family is equipped with the latest Intel Xeon CPUs, so we can try those. Again, configuring this in the Inductiva API is as straightforward as modifying a single argument in the machine's configuration to `c3-standard-88`:

```python
# Configure a dedicated machine from the c3 family boasting 88 vCPUs
my_machine_group = inductiva.resources.MachineGroup(
    machine_type="c3-standard-88")
```

This machine setup further reduces the simulation time to just **5m46s**, achieving
an additionally 1.75x speed-up over the `c2-standard-60` setup. Observe that, in this case, the speed up obtained is higher than the corresponding increase in the number of vCPUs with respect to `c2`-based VMs. This because `c3` VMs are supported by much more recent hardware. And, of course, access to better hardware has a cost.

## Show me the money...

It's important to consider the cost of using these exclusive and more powerful machines from the Google Cloud Platform. At the time of writing, the cost per hour for the `c2-standard-60` is \$3.446, and for the more powerful `c3-standard-88` is \$5.053. Now, compare these values to only \$0.23
per hour for the `c2-standard-4` VMs.

So, while we aim for speed, we must also be mindful of the price of the VMs used, as it directly influences the overall cost of generating our dataset. In the following table, we can see a snapshot of how different machine configurations impact both performance and cost. In the last common, we show the cost of running 1 simulation, that is, the total cost of utilizing a VMs for the time required to run the (same) simulation from beginning to end.

| Machine Type     | Time to run | Hourly rate | Cost of 1 simulation |
| ---------------- | ----------- | ----------- | -------------------- |
| `c2-standard-4`  | 29m27s      | 0.23 \$     | 0.11 \$              |
| `c2-standard-60` | 10m07s      | 3.446 \$    | 0.58 \$              |
| `c3-standard-88` | 5m46s       | 5.053 \$    | 0.49 \$              |

While the `c2-standard-4` machine offers the lowest relative and absolute cost, making it an economical choice for less demanding tasks or if you can simply wait for long enough. On the
other hand, the `c3` machine category stands out for its excellent performance, though
this comes at an absolute cost approximately 4.5 times higher than that of the `c2-standard-4`.

In computation you tend to pay (disproportionately) for speed. Because cost is such an important aspect of creating a large dataset, we will dive deeper on this issue later.

But first, we need to find a way of manipulating the (hyper)parameters of the
simulation programmatically, so we can easily generate thousands of variations of our "base case". This can be done with the templating mechanism provided by the API.
