![Inductiva_Logo_noBG(259 x 76 px)](https://github.com/inductiva/inductiva/assets/102880975/c67e31eb-2775-4ee5-ac64-027f4591c0cc)
# Inductiva: a Python package for running large-scale simulations of physical phenomena on the Cloud.

Inductiva is a Python package designed for researchers, scientists, and engineers 
focused on leveraging the full potential of [free and open-source physical simulation software packages](). 
It enables users to perform complex simulations without local hardware constraints, 
focusing on improving computational efficiency rather than expertise in complex 
systems designed for scaling simulations. 

This introduction is designed to provide you with an overview on Inductiva API, 
you will learn:

* [What is the Inductiva API?]()
* [What are Inductiva API's core functionalities?]()
* [What can you do?]()

## What is the Inductiva API?

Inductiva API is a Python package that enhances the usability and scalability of 
[free and open-source physical simulation software packages](), which are pivotal 
in advancing science and engineering. These software packages, while powerful, often 
present significant challenges in terms of installation, operation, and particularly, 
scaling for large-scale computations. This complexity can restrict the potential 
of engineers and scientists to fully exploit these tools.

Inductiva addresses these challenges by providing on-demand scaling capabilities 
to high-performance systems, making it simpler for users to deploy these simulators 
at scale without the necessity of specialized systems engineering or high-performance 
computing (HPC) team support. 

In essence, **Inductiva's focus lies in computation**: making advanced [computational resources]()
user-friendly and within reach, opening up new possibilities in the field of 
simulation-based research and analysis.

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
simulations in the cloud, leveraging the scalability and power of Inductiva's 
[computational resources](), while bypassing the limitations of your local hardware.

This API offers a range of benefits for diverse simulation needs:

|  	|  	|
|---	|---	|
| Unified Platform 	| _It unifies multiple simulation packages on a single platform, covering various <br>multi-physics aspects that can encompass domains like fluid dynamics, molecular dynamics, <br>plasmas, and structural mechanics. This integration allows for a streamlined, user-friendly <br>access point, eliminating the need to manage separate software for each type of simulation._ 	|
| Hassle-Free Setup 	| _It eliminates the need for installing and managing <br>complex simulation software and corresponding dependencies_ 	|
| Python-Powered Flexibility 	| _It empowers you to write simple <br>Python scripts that blend seamlessly with your existing codebase and machine learning <br>frameworks, unlocking a world of customization and control._ 	|
| Smart Hardware Optimization 	| _It automatically tunes hardware configurations – <br>whether you need CPU or GPU, decides the right number of cores, RAM, etc., specific <br>to each simulation type._ 	|
| Simulate at Scale 	| _It allows running hundreds or even thousands of <br>simulations simultaneously without complex coding from <br>your end._ 	|
 
## What can you do with Inductiva API?

Inductiva's API offers you two distinct options for computational resources:

* [Standard Machines](): Start simulations on open-source simulators quickly using our provided machines with your existing configuration files. These machines handle basic tasks well but aren't suited for high-demand simulations.

* [MachineGroup](): For larger projects needing more power, you can set up a ['MachineGroup'](). This lets you customize a set of virtual machines for heavy-duty simulations, giving you the control and capacity for complex tasks.

----

**Are you ready to scale up your simulations?**
[Getting started](https://inductiva-research-labs-inductiva.readthedocs-hosted.com/en/development/Install.html#) with Inductiva takes only a few minutes.

