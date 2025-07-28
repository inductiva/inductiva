# Versions and Containers 🛠️

## About COAWST
The [COAWST](https://www.usgs.gov/centers/whcmsc/science/coawst-a-coupled-ocean-atmosphere-wave-sediment-transport-modeling-system) COAWST Modeling System is an open-source tool that combines an ocean model, an atmosphere model, a wave model, and a sediment transport model for coastal change studies.

Specifically, COAWST includes an ocean component (ROMS), an atmosphere component (WRF), a hydrology component (WRF_Hydro), wave components (SWAN, WW3 and InWave), a sediment component (USGS Community Sediment Models) and a sea ice model.

## Supported Versions
Inductiva stays up to date with the latest versions of COAWST. Below is a list of the supported versions, along with their respective release dates:

- **3.8** (Jul., 2023)

If you need to use a version not listed here, please feel free to [Contact Us](mailto:support@inductiva.ai).
We’ll be happy to accommodate your request!

## Container Images
Each version of COAWST in the Inductiva API has its own publicly available container image, 
so you can also use it to run simulations. These images are hosted in our Docker Hub repository, 
[Kutu](https://hub.docker.com/r/inductiva/kutu/tags?name=coawst), and you can find the 
Dockerfile details for each version [here](https://github.com/inductiva/kutu/tree/main/simulators/coawst).