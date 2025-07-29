# Versions and Containers 🛠️

## About FVCOM
[FVCOM](https://github.com/FVCOM-GitHub/FVCOM) (Finite Volume Community Ocean Model) is a 3D hydrodynamic model specifically designed for simulating coastal and ocean dynamics. It uses an unstructured grid and finite-volume methods, making it highly adaptable for modeling complex coastlines, estuaries, and bathymetry. It’s a valuable tool for studying marine ecosystems, coastal environments, and the impacts of climate change on ocean systems.

## Supported Versions
Inductiva stays up to date with the latest versions of FVCOM. Below is a list of the supported versions, along with their respective release dates:

- **5.1.0** (Feb., 2023) 

Inductiva has compiled two FVCOM builds:
- **fvcom**: The standard version for general use
- **fvcom_estuary**: Configured to run the Estuary test case included with FVCOM

If you need to use a version not listed here, please feel free to [Contact Us](mailto:support@inductiva.ai).
We’ll be happy to accommodate your request!

## Container Images
Each version of FVCOM in the Inductiva API has its own publicly available container image, 
so you can also use it to run simulations. These images are hosted in our Docker Hub repository, 
[Kutu](https://hub.docker.com/r/inductiva/kutu/tags?name=fvcom), and you can find the 
Dockerfile details for each version [here](https://github.com/inductiva/kutu/tree/main/simulators/fvcom).