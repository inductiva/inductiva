# Introduction

At a very high-level, the functionalities of the API are quite simple. Using simple
python scripting from your laptop, the API allows you:

1. start remote VMs (and shut them down, etc), either working independently or as a
single large MPICluster, with your favorite simulation software already pre-installed;
2. send your simulation scripts from your laptop to one or more remote machines,
and start the simulators;
3. download the results of the simulation, either every data file produced or just
some specific files.

That’s it! Of course, there are many details and possible variations, but this is
the basic flow of usage of the API. So, this section will give you an overview of
how the API works, going into a bit more detail into how the API operates under the
hood. 

More specifically, we will introduce you to:

- [Tasks](./tasks). This is the core computational object of the platform. Your simulation requests are Tasks to be executed remotely on the computational infrastructure we make available
- [Shared and Dedicated Resources](./computational_resources_overview). Tasks get executed on computational resources, which can be Shared by all users, or which can be Dedicated to singe users
- [Storage and Data Flow](./data_flow). The platforms move a lot of data around. Your input files need to be sent to remote computational resources and the output of the simulations needs to be made available to you. 
- [Infrastructure](./infrastructure). Under the hood, there is an infrastructure that provides computational power for executing your tasks. We currently use Google Cloud Platform, and very soon we will have other options available.
