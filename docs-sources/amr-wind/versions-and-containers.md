# ⚙️ Versions and Containers

## About AMR-Wind
[AMR-Wind](https://github.com/Exawind/amr-wind) is a massively parallel, block-structured adaptive-mesh, incompressible flow solver for wind turbine and wind farm simulations. The primary applications for AMR-Wind are: performing large-eddy simulations (LES) of atmospheric boundary layer (ABL) flows and simulating wind farm turbine-wake interactions using actuator disk or actuator line models for turbines.

## Supported Versions
Inductiva stays up to date with the latest versions of AMR-Wind. Below is a list of the supported versions, along with their respective release dates:

- **3.7.0** (Sep., 2025) - with GPU support
- **3.6.0** (Jul., 2025) - with GPU support
- **3.5.0** (Jun., 2025) - with GPU support
- **3.4.1** (Apr., 2025) - with GPU support
- **3.4.0** (Feb., 2025) - with GPU support
- **1.4.0** (Apr., 2024) - with GPU support

If you need to use a version not listed here, please feel free to [Contact Us](mailto:support@inductiva.ai).
We’ll be happy to accommodate your request!

## Container Images
Each version of AMR-Wind in the Inductiva API has its own publicly available container image, 
so you can also use it to run simulations. These images are hosted in our Docker Hub repository, 
[Kutu](https://hub.docker.com/r/inductiva/kutu/tags?name=amr-wind), and you can find the 
Dockerfile details for each version [here](https://github.com/inductiva/kutu/tree/main/simulators/amr-wind).