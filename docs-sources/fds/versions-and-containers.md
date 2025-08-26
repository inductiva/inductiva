# ⚙️ Versions and Containers

## About FDS
[FDS](https://pages.nist.gov/fds-smv/) (Fire Dynamics Simulator) is a computational fluid dynamics (CFD) model designed to simulate fire-driven fluid flows. FDS numerically solves the Navier-Stokes equations to model low-speed, thermally-driven flows with a focus on smoke and heat transport from fires. It is commonly used in simulations for designing smoke management systems, evaluating sprinkler and detector activation, and fire reconstructions.

## Supported Versions
Inductiva stays up to date with the latest versions of FDS. Below is a list of the supported versions, along with their respective release dates:

- **6.10.1** (Apr., 2025) 
- **6.9.1** (Apr., 2024)
- **6.8** (Apr., 2023) 

If you need to use a version not listed here, please feel free to [Contact Us](mailto:support@inductiva.ai).
We’ll be happy to accommodate your request!

## Container Images
Each version of FDS in the Inductiva API has its own publicly available container image, 
so you can also use it to run simulations. These images are hosted in our Docker Hub repository, 
[Kutu](https://hub.docker.com/r/inductiva/kutu/tags?name=fds), and you can find the 
Dockerfile details for each version [here](https://github.com/inductiva/kutu/tree/main/simulators/fds).