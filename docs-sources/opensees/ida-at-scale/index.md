# IDA at Scale: Seismic Assessment of URM Structures

*This tutorial was written by* [Daniel Caicedo](mailto:caicedod93@gmail.com) *in collaboration with the* **Inductiva Team**

<br>

**Incremental Dynamic Analysis (IDA)** is a computational method widely used in **Performance-Based Earthquake Engineering (PBEE)** to assess the seismic performance of structures. It generates response curves across a range of earthquake intensity levels, enabling a probabilistic evaluation of structural behaviour and associated seismic risk.

However, applying this method to **unreinforced masonry (URM)** buildings has historically been limited due to two main challenges:
1. The **complexity of numerical models** needed to accurately capture URM behavior;
2. The **computational burden** of executing hundreds or thousands of simulations.

The first challenge is addressed by developing efficient structural models in OpenSees, using three-
dimensional macroelements to capture the in-plane (IP) and out-of-plane (OOP) mechanisms. These models also account for non-linear floor-to-wall connections and wall-to-wall interlocking.

The second, and often more limiting, challenge lies not only in running a large number of simulations, but also in having the infrastructure required to support them. 

Running simulations on cloud-based HPC is not required to perform IDA, but it offers clear advantages. It removes the need to 
manage physical infrastructure, provides dedicated compute resources, and enables parallel execution at scale. This reduces overall 
runtime and makes large-scale studies more efficient and accessible.

This is where **Inductiva** steps in, allowing you to focus on developing accurate structural models and analysing results, without 
being constrained by **hardware limitations** or **sequential processing**.

In this tutorial, we walk through a practical example involving a two-storey unreinforced masonry building, representative of the 
pre-code masonry building stock in the Lisbon metropolitan area.

All data used in this example was collected as part of the [STAND4HERITAGE](https://stand4heritage.org) project — New STANDards for 
seismic assessment of built cultural HERITAGE.

<p align="center"><img src="../_static/image1" alt="caption1" width="700"></p> - ask Daniel

The tutorial is divided into the following sections:
- [Kickstart Your Simulation Workflow](sections/section1)
- [Run 50 Simulations in Parallel](sections/section2)
- [Results and Key Takeaways](sections/section3)

```{toctree}
:hidden:
sections/section1.md
sections/section2.md
sections/section3.md
```