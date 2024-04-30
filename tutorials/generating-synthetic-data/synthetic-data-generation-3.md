---
myst:
  html_meta:
    description: "Explore available computational resources through the Inductiva API and learn how to choose more powerful hardware to run the 'base case'"
    keywords: "Inductiva API, Programming, HPC, Simulation, Tutorial, Synthetic Data Generation, Physics-ML, SPH"
---

# Choosing your Hardware

## Understanding what is happening under the hood
By default, any task submitted via the API is sent to what we call the “default queue” for execution. The default queue is a shared resource pool of virtual machines (VMs )accessible to all API users. This default queue offers an easy and affordable way for users to prototype and run small-scale simulations similar to our initial low particle count simulation, without the hassle of hacing to manage any computational resource. The goal of the default queue is to provide this feeling of magic, because things will just work out of the box with minimal configuration effort. 

However, because of its shared nature 
and finite capacity, any job sent to the default queue can take a long time to be picked up and executed. This is especially so for higher fidelity simulations rquiring many computation cycles, as it was evident when our enhanced simulation took 30 minutes to run.

To speed up our simulation, our API provides a mechanism that allows us to **launch dedicated machines** exclusively for our projects. These dedicated machines can offer significantly more compute power compared to the standard options accessible through the default queue.

We can choose a setup that aligns best with our simulation needs by instantiating a MachineGroup and sending out simulation task to be executed there. This will make it easier for us to adjust our "base case" parameters, allowing for a higher
particle count that aligns with the more complex scenarios described in the paper, while still having the simulations being executed in a reasoable amount of time.

## Browsing our Available Machines

Our default queue uses virtual machines VMs from Google Cloud Platform, specifically `c2-standard-4.` VMs, based on Intel Xeon Cascade Lake (2nd Gen) processors. These VMs only have 4 vCPUs and, as we saw before, will require about 30 minutes to execute our simulations. The best way to speed things up is to use our API to launch more powerful user-dedicated machines, reserved for our own exclusive use. 

Let's explore alternative machine setups avialable via de APi by running a simple command through the Inductiva [Command Line Interface (CLI)](https://docs.inductiva.ai/en/latest/cli/cli-overview.html):

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

We'll run our "base case" with a _particle radius of 0.008_ again, but now we'll boost our computing power by scaling up the number of vCPUs in the `c2` family available in the default queue. We're moving from the 4 vCPUs of the `c2-standard-4` to a 
much heftier `c2-standard-60` setup, which hopefully will allow us to significantly decrease the simulation time. This can easily be achieved with just a few extra lines of code:

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
task.download_outputs()

# Terminate the machine group
my_machine_group.terminate()
```
Notice the significant reduction in runtime for this simulation: down from 30m 
to just **10m07s**, achieving an approximately 3-fold increase in speed! Remember that, in practice, speed ups are not linear with the number of vCPUs, so this is alreasy quite an achievement with just a few lines of code.

We can enhance the performance further by choosing even more powerful machines. For example, the `c3` family is equipped with the latest Intel Xeon CPUs, so we can try those. Again, configuring this 
in the Inductiva API is as straightforward as modifying a single argument in  
the machine's configuration to `c3-standard-88`:

```python
# Configure a dedicated machine from the c3 family boasting 88 vCPUs
my_machine_group = inductiva.resources.MachineGroup(
    machine_type="c3-standard-88")
```
This machine setup further reduces the simulation time to just **5m46s**, achieving 
an aditionally 1.75x speed-up over the `c2-standard-60` setup. Observe that, in this case, the speed up obtained is higher than the corresponding increase in the number of vCPUs with respect to `c2`-based VMs because `c3` VMs are supported by much more recent hardware.

## Show me the money...
Of course, there’s a catch: these exclusive machines have a certain non-negligible cost per
hour. It's important to consider the cost of using these exclusive, more powerful 
machines. At the time of writing, the cost per hour for the `c2-standard-60` is $3.446, 
and for the more powerful `c3-standard-88`, it's $5.053, compared to only $0.23 
per hour for the `c2-standard-4` used in the default queue. While we aim for 
efficiency and speed, we must also be mindful of the price, as it directly influences 
the overall cost of generating our dataset.

In the following table, we can see a snapshot of how different machine configurations impact 
both performance and cost.

| Machine Type | Time to run | Hourly rate | Total Cost |
| --- | --- | --- | --- |
| `c2-standard-4` | 29m27s | 0.23 \$ | 0.11 \$ |
| `c2-standard-60` | 10m07s | 3.446 \$ | 0.58 \$ |
| `c3-standard-88` | 5m46s | 5.053 \$ | 0.49 \$ | 

While the `c2-standard-4` machine offers the lowest cost, it takes much longer 
to compute, making it an economical choice for less demanding tasks. On the 
other hand, the `c3` machine category stands out for its fast performance, though 
it comes at a price approximately four times higher than that of the `c2-standard-4`. 
This preliminary comparison offers a peak into a more detailed analysis in 
upcoming tutorials, where we will delve further into optimizing the trade-off 
between simulation speed and operational costs.

## Up Next: Generalizing our Simulation Script with Inductiva’s Templating Engine

In this step, we touched on the importance of choosing the right machine 
setup and how cost considerations play a crucial role. However, this can become 
quite challenging when we start manipulating some of the hyperparameters of the 
simulator without manually tweaking them to find the best dataset generation method.

In the [next phase]({% post_url 2024-03-24-api-synthetic-data-generation-4 %}) of
this tutorial, we will transform our "base case" simulation configuration files into 
a "programmable script" by using Inductiva's [Templating Engine](https://docs.inductiva.ai/en/latest/explore_api/templating.html). 
This script will enable us to programmatically simulate all sorts of variations 
of the base case, each with unique parameter and hyperparameter settings, all through 
Python scripting!
