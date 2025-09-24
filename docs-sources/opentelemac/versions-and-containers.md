# ⚙️ Versions and Containers

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

## Supported Versions
Inductiva stays up to date with the latest versions of OpenTelemac. Below is a list of the supported versions, along with their respective release dates:

- **9.0.0** (Dec., 2024)
- **8p4r1** (Nov., 2023)
- **8p4r0** (Dec., 2022)
- **8p1r2** (Dec., 2021)

If you need to use a version not listed here, please feel free to [Contact Us](mailto:support@inductiva.ai).
We’ll be happy to accommodate your request!

## Container Images
Each version of OpenTelemac in the Inductiva API has its own publicly available container image, 
so you can also use it to run simulations. These images are hosted in our Docker Hub repository, 
[Kutu](https://hub.docker.com/r/inductiva/kutu/tags?name=opentelemac), and you can find the 
Dockerfile details for each version [here](https://github.com/inductiva/kutu/tree/main/simulators/opentelemac).