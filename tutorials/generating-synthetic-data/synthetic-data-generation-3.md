---
myst:
  html_meta:
    description: "Explore available computational resources through the Inductiva API and learn how to choose more powerful hardware to run the 'base case'"
    keywords: "Inductiva API, Programming, HPC, Simulation, Tutorial, Synthetic Data Generation, Physics-ML, SPH"
---

# Choose the Hardware

In our [previous step](synthetic-data-generation-2), 
we learned how to set up a "base case" simulation, building upon the study by [Sanchez-Gonzalez et al.](https://arxiv.org/abs/2002.09405) where a Graph Neural Network (GNN) is used to learn how to simulate fluid dynamics.
Our "base case" simulated a 0.5m cube of water released from one of the top corners 
of a sealed 1m cubic box splashing against the box's walls over a 
4-second duration.

After our initial run of the "base case" simulation with a _0.01_ 
particle radius and later adjusting the particle count to _0.008_ to align with 
the research we're basing this on, we noticed a much longer runtime — specifically, 
**30 minutes for a single run simulating 4 seconds of fluid dynamics**. This means 
that if we want to run 10,000 or even 20,000 different variations to generate a large 
enough dataset, we'd be facing thousands of hours in computation time, and utimately, a need for more 
powerful hardware. This time around, we're looking at how to speed things up with 
beefier hardware through our API, and what that means in terms of costs. Stick around 
as we dive into choosing the right machine for the job and weigh out the price points.

## Choosing more Powerful Hardware

When we submitted our "base case" simulation task, it was sent to what we call the “default queue”,
a queue of tasks processed by a shared resource pool of virtual machines accessible to all 
API users. This default queue offers an easy and affordable way for users to 
prototype and run small-scale simulations similar to our initial low particle count 
simulation, without the hassle of managing resources. However, its shared nature 
and finite capacity mean that simulations, especially higher fidelity simulations, 
can take a long time to be picked up and to execute. This was evident when our 
enhanced simulation took 30 minutes to run on the default queue's virtual 
machines (VMs), the slowest option available, due to their limited computational power.

To speed up our simulation, our API provides a mechanism that allows us to **launch dedicated machines** 
exclusively for our projects. These dedicated machines can offer 
significantly more compute power compared to the standard options accessible through the 
default queue. By instantiating a MachineGroup, we now have the flexibility to 
choose the setup that aligns best with our simulation needs. This will make it
easier for us to adjust our "base case" parameters, enhancing its resolution and 
increasing the particle count to align with the more complex scenarios described 
in the research we're building upon.

### Browsing our Available Machines

Our default queue uses virtual machines from the `c2` family, specifically `c2-standard-4.` 
These only have 4 vCPUs and will take a long time to compute complex simulations.
To speed things up, the best way is to use our API to launch dedicated machines
reserved for our exclusive use. Let's explore our options for machine setups by
running a simple command through the Inductiva [Command Line Interface (CLI)](https://docs.inductiva.ai/en/latest/cli/cli-overview.html) tool:

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
### Running our "Base Case" on Dedicated Machines

We'll run our "base case" with a _particle radius of 0.008_ again, but now we'll 
boost our computing power by scaling up the number of vCPUs in the `c2` family available 
in the default queue. We're moving from the 4 vCPUs of the `c2-standard-4` to a 
much heftier `c2-standard-60` setup, allowing us to significantly increase processing 
power. This can easily be achieved with just a few extra lines of code:

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
to just **10m07s**, achieving an approximately 3-fold increase in speed!

Of course, we can enhance the performance further with more powerful machines equipped 
with the latest Intel Xeon CPUs, specifically the `c3` family. Configuring this 
in the Inductiva API is as straightforward as modifying a single argument in  
the machine's configuration:

```python
# Configure a dedicated machine from the c3 family boasting 88 vCPUs
my_machine_group = inductiva.resources.MachineGroup(
    machine_type="c3-standard-88")
```
This machine setup further reduces the simulation time to just **5m46s**, achieving 
a significant 1.75x speed-up over the `c2-standard-60` setup.

Of course, there’s a catch: these exclusive machines have a non-negligible price per
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
