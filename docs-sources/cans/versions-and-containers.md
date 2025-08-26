# ⚙️ Versions and Containers

## About CaNS
[CaNS](https://github.com/CaNS-World/CaNS) (Canonical Navier-Stokes) is a high-performance simulator for massively-parallel numerical simulations of fluid flows. It is specifically designed to solve problems involving incompressible, Newtonian fluids, making it a popular choice for researchers and engineers in fluid dynamics. CaNS is widely used in simulations of turbulence, vortex dynamics, and large-scale atmospheric or oceanic flows.

## Supported Versions
Inductiva stays up to date with the latest versions of CaNS. Below is a list of the supported versions, along with their respective release dates:

- **3.0.0** (Apr., 2025) - with GPU support
- **2.4.0** (Jan., 2025) - with GPU support
- **2.3.4** (Apr., 2024) for CPU

If you need to use a version not listed here, please feel free to [Contact Us](mailto:support@inductiva.ai).
We’ll be happy to accommodate your request!

## Container Images
Each version of CaNS in the Inductiva API has its own publicly available container image, 
so you can also use it to run simulations. These images are hosted in our Docker Hub repository, 
[Kutu](https://hub.docker.com/r/inductiva/kutu/tags?name=cans), and you can find the 
Dockerfile details for each version [here](https://github.com/inductiva/kutu/tree/main/simulators/cans).