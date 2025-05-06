# Welcome to OpenTelemac at Inductiva
Your resource hub for all things OpenTelemac at Inductiva. Whether you're just starting out or an experienced user, you'll find the resources you need to seamlessly run your OpenTelemac simulations on Cloud machines equipped with hundreds of cores and terabytes of disk space.

Inductiva simplifies research by making high-performance computing more accessible and cost-effective. Use the power of the Cloud to **scale your simulations** and **finish your projects sooner**, while keeping your costs in check! 

## About OpenTelemac
[OpenTelemac](https://www.opentelemac.org) is a suite of open-source numerical models designed to 
simulate free surface water flow, sediment transport, waves and water quality in rivers, 
estuaries and coastal areas. It provides 2D, 3D and multiphase models designed to address 
complex hydrodynamic and environmental challenges. OpenTelemac is widely used for water management, 
flood forecasting, and coastal engineering applications.

OpenTelemac has a modular structure and supports unstructured meshes, high performance computing (HPC) and robust numerical solvers. Its architecture allows the integration of multiple physical processes under a unified computational framework. The suite includes:
- **TELEMAC-2D**: Simulates 2D shallow water flows using the Saint-Venant equations. Ideal for tides, flood modeling, and river hydraulics.
- **TELEMAC-3D**: A 3D non-hydrostatic model solving the Navier-Stokes equations with support for turbulence, stratification and tracer transport.
- **GAIA**: Models sediment transport in 2D and 3D environments, including bedload and suspended sediment processes.
- **TOMAWAC**: Spectral wave model to simulate wave propagation from offshore to coastal areas.
- **ARTEMIS**: Treats wave agitation in harbors and near structures based on mild slope equations.
- **MASCARET**: One-dimensional hydraulic solver ideal for river and canal flows.
- **WAQTEL**: Simulates water quality and thermal processes, optionally coupled with the AED2 model for advanced biochemical simulations.
- **KHIONE**: Focuses on ice formation and water-air-ice interactions in cold regions.

## What You'll Find Here
- **Tutorials:** Step-by-step guides to help you learn how to run OpenTelemac through the Inductiva API. From getting started to advanced tutorials, we have you covered.
- **Benchmarks:** A trusted guide to selecting the right simulation hardware for your needs. These benchmarks, conducted using the Inductiva platform, provide insight into how OpenTelemac performs on different hardware configurations.


```{toctree}
---
caption: " "
maxdepth: 2
hidden: true
---
versions-and-containers
configuring-telemac-simulations
```

```{toctree}
---
caption: Tutorials
maxdepth: 2
hidden: true
---
setup-test
quick-start
```