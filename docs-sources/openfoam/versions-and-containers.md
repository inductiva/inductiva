# Versions and Containers 🛠️

## About OpenFOAM
OpenFOAM is an open-source C++ toolbox designed for the development of custom numerical solvers and pre- and post-processing utilities. It is widely used for solving a variety of continuum mechanics problems, including computational fluid dynamics, solid mechanics, electromagnetism, and multiphase flow.

## Supported Versioms
Inductiva stays up to date with the latest versions of OpenFOAM and supports both OpenFOAM distributions: [ESI Group](https://www.openfoam.com/about-esi-opencfd) and [OpenFOAM Foundation](https://openfoam.org/), giving you access to the full range of OpenFOAM capabilities.

Below is a list of the supported versions, along with their respective release dates:

OpenFOAM (ESI)
- **24.12** (Dec., 2024)
- **24.06** (Jun., 2024)
- **22.06** (Jun., 2022)

OpenFOAM (Foundation)
- **12** (Jul., 2024)
- **8** (Jul., 2020)

If you need to use a version not listed here, please feel free to [Contact Us](mailto:support@inductiva.ai).
We’ll be happy to accommodate your request!

## Container Images
Each version of OpenFOAM in the Inductiva API has its own publicly available container image, 
so you can also use it to run simulations. These images are hosted in our Docker Hub repository, 
[Kutu](https://hub.docker.com/r/inductiva/kutu/tags?name=openfoam), and you can find the 
Dockerfile details for each version [here](https://github.com/inductiva/kutu/tree/main/simulators/openfoam).