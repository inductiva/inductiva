# The Inductiva Guide to REEF3D 🌊
Your resource hub for all things REEF3D at Inductiva. Whether you're just starting out or an experienced user, you'll find the resources you need to seamlessly run your REEF3D simulations on Cloud machines equipped with hundreds of cores and terabytes of disk space.

Inductiva simplifies research by making high-performance computing more accessible and cost-effective. Use the power of the Cloud to **scale your simulations** and **finish your projects sooner**, while keeping your costs in check! 

## About REEF3D
[REEF3D](https://reef3d.wordpress.com/) is an open-source hydrodynamics framework specifically designed for coastal, marine, and hydraulic engineering applications. 
Built with a modular programming approach, it offers multiphysics solvers tailored to address fluid flow problems, such as sediment transport and floating body dynamics, as well as wave modeling.

The modular programming approach allows the framework to incorporate a range of different flow solvers which together represent all relevant length scales. Depending on the wave or flow conditions, the following optimized hydrodynamic modules are available:
- **REEF3D::CFD** solves the Navier-Stokes equations in three dimensions. For near-field simulations with a complex free surface pattern, it uses a two-phase flow approach with the level set method for interface capturing.
- **REEF3D::FNPF** is a three-dimensional fully nonlinear potential flow solver. It is massively parallelized and can be used to create large-scale phase-resolved sea states at all water depths.
- **REEF3D::SFLOW** is a depth-averaged model, solving the non-hydrostatic shallow water equations ideal for near-shore hydrodynamics and river flow.

## What You'll Find Here
- **Tutorials:** Step-by-step guides to help you learn how to run REEF3D through the Inductiva API. From getting started to advanced tutorials, we have you covered.
- **Benchmarks:** A trusted guide to selecting the right simulation hardware for your needs. These benchmarks, conducted using the Inductiva platform, provide insight into how REEF3D performs on different hardware configurations.

```{toctree}
---
caption: " "
maxdepth: 1
hidden: true
---
versions-and-containers
```

```{toctree}
---
caption: Tutorials
maxdepth: 3
hidden: true
---
setup-test
quick-start
run-3d-dam-break-scenario
```


