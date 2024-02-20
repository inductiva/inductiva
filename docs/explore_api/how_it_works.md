# How it works

At a very high-level, the functionalities of the API are quite simple. Using simple 
python scripts from your laptop, you can streamline your overall simulation workflow.

In essence, the API enables you to:

<div align="center">
   <img src="../_static/infographic-apifunctionality-fullscreen.svg" alt="Inductiva API Usage Flow">
      <figcaption align = "center"><b>Inductiva API Usage Flow</b></figcaption>

</div>


1. Start and manage remote Virtual Machines (VMs), operating either independently 
or collectively as an MPICluster, that are equipped with pre-installed simulation 
software.

2. Send your simulation scripts from your laptop to one or more remote machines, 
and start the simulators.

3. Download simulation results, selecting either all data files generated or 
specific ones of interest.

That’s it! While this outlines the basic flow of usage of the API, there are
nuances and additional options available, which we'll explore in detail throughout
this section. More specifically, you'll become familiarized with how the API works 
as we introduce you to:

- [Tasks](./tasks): The core computational object of the platform, representing your 
simulation requests to be processed remotely on the computational infrastructure 
we make available.

- [Shared and Dedicated Resources](./shared_dedicated_resources): Understanding 
how tasks get executed on computational resources that are either **shared** among 
all users or **dedicated** to individual users.

- [Storage and Data Flow](./data_flow): Understanding the flow of data across platforms, from 
sending your input files to remote computational resources to accessing the 
simulation outputs.

- [Infrastructure](./infrastructure): A further look under the hood, representing 
the underlying infrastructure powering your simulations, currently leveraging the 
Google Cloud Platform, with more options on the horizon.

