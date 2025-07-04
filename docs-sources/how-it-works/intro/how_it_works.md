# How It Works

At a very high-level, the functionalities of the API are quite simple. Using simple
python scripts from your laptop, you can streamline your entire simulation workflow.

In essence, the API enables you to:

<div align="center">
   <img src="../_static/infographic-apifunctionality-fullscreen.svg" alt="Inductiva API Usage Flow">
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
this section. More specifically, you'll become familiarized with how the API
works:

- [Tasks](./tasks.md): Learn about tasks, the API's core computational object, which
gets created when you submit your simulation request, enabling you to track its
progress and access its outputs in real-time.

- [Resource Allocation Options](./shared_dedicated_resources.md): Discover how
the API allocates and runs your simulation tasks on **dedicated** computational resources.

- [Storage and Data Flow](./data_flow.md): Get to know the typical flow of data
when you invoke a remote simulator using the Inductiva API, from sending your input
files to remote computational resources to accessing the simulation outputs.

- [Computational Infrastructure](./computational-infrastructure.md): Take a further
look under the hood to learn the underlying infrastructure powering your simulations,
especially how the API enables you to access a variety of computational resources,
manages your computational workload, and allocates it to the appropriate computational
resource through a unified Python code.

- [Benchmarking](./benchmarking.md): Evaluate and compare the performance and cost
of different machine configurations for your computational jobs by running a sample
simulation across a variety of Virtual Machine (VM) options and get actionable insights.

- [Templating Engine](./templating.md): Explore how the API enables you to transform
fixed parameter values in your "base case" simulation configuration files into
programmable variables you can adjust for running large variations of your simulation.

- [Configuring Simulators](./configuring-simulators.md): Explore how the API
wraps existing open-source software packages within layers that facilitate their
execution across various cloud-based virtual machines and providers, transforming
them into abstract computational tasks.
