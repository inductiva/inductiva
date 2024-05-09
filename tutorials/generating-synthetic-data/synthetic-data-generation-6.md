---
myst:
  html_meta:
    keywords: "Inductiva API, Programming, HPC, Simulation"
orphan: true
---

# Benchmark Computational Resources

In the previous step of this tutorial, we used Inductiva’s templating mechanism to 
transform the configuration files for our base case simulation into a generalized 
file whose parameters we could programmatically change. This way, we can systematically 
produce many variations of the base case we are interested in so that we can then 
simulate to generate a diverse enough dataset for training our ML models. The 
parameters that we transformed into variables were the dimensions, initial position and 
initial velocity of the fluid block, as well as the density and viscosity of the fluid 
itself. These parameters control the core physical aspects of the simulation, and we 
will need to cover a wide range of values over these parameters to create a rich 
dataset. 

Additionally, we also generalized a few hyperparameters that control how the simulation is performed and the level of fidelity, namely: ***particle radius***.

- How much simulation time do we want to generate (e.g. 5 seconds of fluid splashing)?
- How small are the particles that we use to represent the fluid?        

These hyperparameters are not related to the physical aspects of the simulation and
only affect the underlying algorithm to improve the numerical stability and resolution.
Hence, they have a deep impact on how much computing power is required for running a 
simulation. They also affect the amount of data that each simulation run produces. So, 
from a practical standpoint, and assuming that we have a fixed budget for investing in 
the generation of data (e.g. 1000\$ to cover the compute costs), the way we set these 
hyperparameters determines how many simulations we will be able to run and 
how long we will wait until the desired data becomes available. 

Remember that we had concluded that the choice of the machines has an impact on 
how fast the simulation would run, and how much that simulation would cost. So, besides 
the need to assess the impact of the hyperparameters in the computational cost of the 
simulation, we will also need to see if we can optimize the machine type. That’s a lot 
of variables to play with. But first, let’s study the impact of the hyperparameters on 
the running time of a simulation.

# Assessing the impact of changing simulation hyperparameters

Here, we will study varying the number of each of the hyperparameters separately, to see how they impact the total running time of the simulation. Of the 2 hyperparameters, the size of the particles is the one that should affect the simulation time the most: the smaller the particle size, the larger the number of particles. Observe that, to fill the same volume, the number of particles grows cubically as their size shrinks. This also means that the amount of information to process and update grows cubically as the particle size is reduced. Let’s start by testing that, and we can do it using the Inductiva API with a simple for loop. We will keep the amount of simulation time fixed (4 secs). Here is the code:

And these are the results:

[Sérgio acho que podemos pegar numa máquina boa, high-cpu com um bom número de cores e correr para uns 10 ou 15 valores de particle size, se calhar em escala logarítmica… 
No paper eles usaram simulacoes com varios tamanhos, mas a maior tinha 20k, por isso acho que devemos escolher um particle size que nao nos deixe ir muito acima. Na verdade, se for muito maior depois o utilizador final teria problemas em usar GPUs para treinar GNN.

Acho que tb basta 5 segundos de tempo simulado] [temos info sobre o número de partículas usadas em cada simulação? se sim era fixe tb ter isso no gráfico. Não sei se é propriedade que fica acessível só nos logs. Valeria a pena que isso ficasse disponível num (sub)campo do objecto task?]
 
So, as we can see, the computation times grow pretty quickly. And the prices for each simulation go up accordingly. If we run simulations with a particle size of XX (quite small, high fidelity), we should expect them to cost XX using the chosen hardware. Therefore, using our budget, we will only be able to run XXX simulations. This is not a lot if we think about the number of parameters relevant to the physical properties of the use case that we can play with. As we will see later, there is a smart way of saving costs, but let’s proceed with studying the other hyperparameters.

Using very similar code, we can also study the the impact of changing the time step from 0.001 to 0.05 (in logarithmic scale):


And finally, here is the impact of changing the total simulated time, linearly, from 1 to 10 seconds:


So, clearly, the most significant hyperparameter at play is the particle size. Let’s start by finding a good compromise on that one. First let’s remember that in the work by Sanchez-Gonzalez et al. [https://arxiv.org/abs/2002.09405], the simulation used as training data with the largest number of particles had 20k particles (see Appendix B.1). But most simulations were actually much smaller, in the range 1k to 8k. So, if we want to limit ourselves approximately to that range, and even considering that the size of the fluid block can be bigger that the one we are testing now (remember the user the current dimensions), we can keep particle size at XX.

Also, in the same work, the time step chosen was 2.5 ms, generating 400 data frames per second of simulation. We see that, on the time step scale of 1 ms to 5 ms, the simulation time grows in ???. (Should we go to 1ms, or maybe this will generate too much data?). 

Also, Sanchez-Gonzalez et al., simulations have between 1 and 5 seconds of simulated time.  To give users the option of using a little bit more information to train their models on, we will run all the simulation for 5 seconds.

So, choosing these hyperparameters, our base case simulation will run in XX using a CCCC machine. So, assuming variations of this would take approximately this time (some would run in less or more depending on the size of the fluid block), then we would be able to simulate XXX cases on this hardware option. This is still not that many, so we need to go deeper on our study.

# Choosing the best VM for the job

Now that we have found a reasonable instantiation for the simulation hyperparameters, we can try to optimize for the type of VM that we should use. Inductiva’s computation is currently on Google Cloud VMs, and there are dozens of options available at different costs. GCP makes available VMs from multiple family types, e.g. c2, c2d, c3, etc . Different VM families run on different physical hardware, with different CPUs, RAM specs, motherboards, etc. 

Some simulation software was written to benefit from running on some specific hardware types, so their performance may change significantly from one VM family to another. In some cases, running simulations on less expensive VMs may actually lead to better performance than using more expensive VMs. Also, VMs come with a configurable number of vCPUs, which can go from mere 4 vCPUs to several dozens or hundreds as is the case of c2d and c3d machines. Typically, the more vCPUs you use the faster your simulation runs, but obviously, since scaling is never linear, there is a sweet spot. Again, using a machine with a a larger number of vCPU may not pay off since the efficiency of the parallelization decreases and, so the price you pay rises faster than the speed up benefits you collect.

So, if we are running thousands of simulations, we better find the most cost-effective hardware so that we can run as many simulations as possible within our limited budget. Fortunately, Inductica API allows you to easily obtain performance statistics over all hardware options available. In short, the API allows you to get the list of all available VM options, and you can simply iterate through all of them using a “for loop” and collect the statistics you need. This is how you would run the base case, with the previously chosen hyperparameters, over the various vCPU configurations for 5 VMs families:


With a bit of post-processing, you will get the following performance chart:

[time to execute 1 base case simulation / vCPu for each family]

And since we care mostly about cost (assuming we are ok with waiting a some time for the simulations to run) here this performance mapped to cost of simulation:

[cost of 1 base case simulation / vCPU for each family]

Please note that these costs have been computed using the prices made available by GCP for each VM type at the time of writing. So, we conclude that we would be better off using VMs of one of the following types:

A  - some description and approximate cost per simulation
B
C

These would allow us to run about ZZ simulations using the XXX budget that we have available for data generation. But there is still one more lever we can pull: using “spot instances” instead of “dedicated instances”. 

# Spot Instances: one more cost saving strategy
Spot machines are VM that are available for computations as long as no one else requests them as dedicated instances. They can be seen as spare capacity that is not being requested by any user, so the cloud provider is better off renting them at a much cheaper price to anyone who is ok with losing access to them at any moment. So, running loads on spot instances does not guarantee that your job will get to the end. So, you may spend money executing part of the computation and then lose everything close to the end. 

However, the price at which you get those machines is 3 to 5 times lower than that of dedicated machines, which means that if you are ok with losing jobs once in a while, and that happens rarely (e.g. less than 1 in 4 times), then it is much cheaper to run them in a spot instance. This is especially so, if you can restart your jobs and try again automatically, as happens if you submit the job using the Inductiva API. 

Below, we show how you can request your job to run on a spot instance instead of running on a dedicated instance. You basically only need to set the parameter spot=True when starting your Machine Group:


Thus, we can assume that the total cost per simulation can be only XXX (again, according to GCP pricing at the time of writing). This will allow us to run approximately XXX simulations, which seems already enough to generate a robust enough training set for the GNN models, especially given the fact that the authors of the paper used about 22k runs for all their tests.

In the next post, we will explain how to actually run all these XXX simulations with only a few lines of code, while distributing the load over XXX VMs in parallel. 


## Varying the particle radius


| Particle Radius | Number of Particles | Size of the data produces / MB |
|-----------------|---------------------|--------------------------------|
| 0.008           | 29791               | 316                            |
| 0.0082          | 27000               | 287                            |
| 0.0085          | 24389               | 259                            |
| 0.0088          | 21952               | 233                            |
| 0.009           | 21948               | 233                            |
| 0.0094          | 19601               | 208                            |
| 0.0097          | 17572               | 186                            |

c3-standard-4 machine

| Particle Radius | Running time              | Cost / USD |
|-----------------|---------------------------|------------|
| 0.008           | 21 minutes 33 seconds     | 0.0825     |
| 0.0082          | 22 minutes 4 seconds      | 0.0845     |
| 0.0085          | 16 minutes 39 seconds     | 0.0638     |
| 0.0088          | 16 minutes 20 seconds     | 0.0625     |
| 0.009           | 17 minutes 46 seconds     | 0.0680     |
| 0.0094          | 13 minutes 15 seconds     | 0.0507     |
| 0.0097          | 13 minutes 32 seconds     | 0.0519     |
| **0.01**        | **11 minutes 39 seconds** | **0.0446** |

c3-standard-8 machine

| Particle Radius | Running time              | Cost / USD |
|-----------------|---------------------------|------------|
| 0.008           | 14 minutes 11 seconds     | 0.1086     |
| 0.0082          | 13 minutes                | 0.0995     |
| 0.0085          | 12 minutes 55 seconds     | 0.0990     |
| 0.0088          | 11 minutes 8 seconds      | 0.0853     |
| 0.009           | 9 minutes 28 seconds      | 0.0726     |
| 0.0094          | 9 minutes 12 seconds      | 0.0705     |
| 0.0097          | 8 minutes 7 seconds       | 0.0623     |
| **0.01**        | ***7 minutes 12 seconds** | **0.0552** |

c3-standard-88

| Particle Radius | Running time            | Cost / USD |
|-----------------|-------------------------|------------|
| 0.008           | 4 minutes 41 seconds    | 0.3957     |
| 0.0082          | 5 minutes 2 seconds     | 0.4244     |
| 0.0085          | 3 minutes 50 seconds    | 0.3234     |
| 0.0088          | 4 minutes 2 seconds     | 0.3407     |
| 0.009           | 3 minutes 37 seconds    | 0.3051     |
| 0.0094          | 3 minutes 42 seconds    | 0.3122     |
| 0.0097          | 3 minutes 19 seconds    | 0.2804     |
| **0.01**        | **3 minutes 4 seconds** | **0.2593** |



## Varying simulation time

**Maybe Change** *running time* to *computation time*

On a c3-standard-4 machine

| Simulation Time / s | Running Time              | Cost       |
|---------------------|---------------------------|------------|
| 1                   | 2 minutes 43 seconds      | 0.0104     |
| 2                   | 5 minutes                 | 0.0191     |
| 3                   | 8 minutes 8 seconds       | 0.0311     |
| **4**               | **11 minutes 39 seconds** | **0.0446** |
| 5                   | 14 minutes 29 seconds     | 0.0555     |

On a c3-standard-8

| Simulation Time / s | Running Time              | Cost       |
|---------------------|---------------------------|------------|
| 1                   | 1 minute 34 seconds       | 0.0121     |
| 2                   | 3 minutes 24 seconds      | 0.0261     |
| 3                   | 4 minutes 41 seconds      | 0.0359     |
| **4**               | ***7 minutes 12 seconds** | **0.0552** |
| 5                   | 8 minutes 36 seconds      | 0.0658     |

On a c3-standard-88

| Simulation Time / s | Running Time            | Cost       |
|---------------------|-------------------------|------------|
| 1                   | 0 minutes 46 seconds    | 0.0648     |
| 2                   | 1 minute 21 seconds     | 0.1144     |
| 3                   | 2 minutes 6 seconds     | 0.1773     |
| **4**               | **3 minutes 4 seconds** | **0.2593** |
| 5                   | 3 minutes 28 seconds    | 0.2932     |


## Varying time step interval

For a c3-standard-4

| Simulation time step interval | Running Time              | Cost       |
|-------------------------------|---------------------------|------------|
| 0.001                         | 16 minutes 22 seconds     | 0.0627     |
| **0.01**                      | **11 minutes 39 seconds** | **0.0446** |
| 0.0255                        | 12 minutes 29 seconds     | 0.0478     |
| 0.05                          | 18 minutes 19 seconds     | 0.0702     |

For a c3-standard-8

| Simulation time step interval | Running Time              | Cost       |
|-------------------------------|---------------------------|------------|
| 0.001                         | 8 minutes 37 seconds      | 0.0660     |
| **0.01**                      | ***7 minutes 12 seconds** | **0.0552** |
| 0.0255                        | 6 minutes 49 seconds      | 0.0522     |
| 0.05                          | 11 minutes 19 seconds     | 0.0867     |

For a c3-standard-88

| Simulation time step interval | Running Time            | Cost       |
|-------------------------------|-------------------------|------------|
| 0.001                         | 4 minutes 28 seconds    | 0.3768     |
| **0.01**                      | **3 minutes 4 seconds** | **0.2593** |
| 0.0255                        | 3 minutes 20 seconds    | 0.2808     |
| 0.05                          | 8 minutes 17 seconds    | 0.6983     |

**Not sure why the values for 0.05 times step are larger.** Maybe the simulation is erroneous and all particles clump together.
