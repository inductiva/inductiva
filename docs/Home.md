![Inductiva_Logo_noBG(259 x 76 px)](https://github.com/inductiva/inductiva/assets/102880975/c67e31eb-2775-4ee5-ac64-027f4591c0cc)
# Inductiva: a Python package for scaling simulations on the Cloud or On-Premise

Inductiva is an API that empowers scientists and engineers to maximize the 
capabilities of [free and open-source simulation software](), enabling the scaling 
of simulation projects across hundreds of machines with just a simple Python script. 
It simplifies the complexities of cloud resource management, cost optimization, and 
configuration, offering a straightforward Python interface for deploying computational 
resources, including MPI clusters. Inductiva streamlines the orchestration of 
simulation pipelines, making advanced simulations more accessible and efficient.

This introduction provides you with an overview on Inductiva API, 
you will get to know:

* [What is the Inductiva API?]()
* [What are Inductiva API's core functionalities?]()
* [What can you do?]()

## What is the Inductiva API?

Inductiva API is a Python package that enhances the usability and scalability of 
[free and open-source physical simulation software packages](). These software 
packages, while powerful, often 
present significant challenges in terms of installation, operation, and particularly 
scaling to large computational infrastructures, such as the Cloud or On-Premise 
clusters. This complexity can restrict the potential 
of engineers and scientists to fully exploit these tools.

Inductiva facilitates this process by acting as an intermediary between users and high-performance computational resources. It simplifies the task of scaling, eliminating the need for specialized systems engineering or an in-house high-performance computing (HPC) team. 

In essence, **Inductiva's focus lies in simplifying the underlying computational infrastructure**, 
allowing users to focus on their simulations without worrying about the complexities 
of computational resources, opening up new possibilities in the field of simulation-based 
research and analysis.

## What are Inductiva API's core funtionalities?
Inductiva designed the core funtionalities of its API tool to facilitate 
various tasks in a computational environment with minimal coding effort. It will 
allow you to:
* easily dispatch 
computational tasks,
* manage system resources,
* and run detailed simulations.

You will interact with these cloud-based simulations through a Python client on 
your local computer. This interface enables you to select, manage, and scale your 
simulations on the Cloud or On-Premise.

This API offers a range of benefits for diverse simulation needs:

|  	|  	|
|---	|---	|
| **Unified Platform** 	| _It unifies standard open-source simulation packages on a single platform, covering various multi-physics aspects that can encompass domains like fluid dynamics, molecular dynamics, plasmas, and structural mechanics. This integration allows for a streamlined, user-friendly access point, eliminating the need to manage separate software for each type of simulation._ 	|
| **Hassle-Free Setup** 	| _It eliminates the need for installing and managing complex simulation software and corresponding dependencies_ 	|
| **Python-Powered Flexibility** 	| _It empowers you to write simple Python scripts that blend seamlessly with your existing codebase and machine learning frameworks_ 	|
| **Smart Hardware Optimization** 	| _It automatically tunes hardware configurations – whether you need CPU or GPU, decides the right number of cores, RAM, etc., specific to each simulation type._ 	|
| **Efficient Simulation at Scale** 	| _It simplifies orchestrating and managing simulations, from auto-scaling to automatic termination, without complex coding from your end. It offers essential features like MPI cluster setup and cost visibility, making it easy to run extensive simulations while controlling expenses._ 	|
 
## What can you do with Inductiva API?

Inductiva's API offers you two distinct options for computational resources:

* [Standard Machines](): Start simulations on open-source simulators quickly using 
our provided machines with your existing configuration files. These machines handle 
basic tasks well but aren't suited for high-demand simulations.

* [MachineGroup](): For larger projects needing more power, you can set up a ['MachineGroup'](). 
This lets you customize a set of virtual machines for heavy-duty simulations, giving 
you the control and capacity for complex tasks.

----

**Are you ready to scale up your simulations?**
[Getting started](https://inductiva-research-labs-inductiva.readthedocs-hosted.com/en/development/Install.html#) with Inductiva takes only a few minutes.

