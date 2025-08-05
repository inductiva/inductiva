# ⚙️ Versions and Containers

## About OpenSees
[OpenSees](https://opensees.berkeley.edu) (Open System for Earthquake Engineering Simulation) is a software framework designed for the development of sequential, 
parallel and grid-enabled finite element applications in earthquake engineering. It allows users to simulate the response of structural and 
geotechnical systems subjected to earthquakes and other hazards using scripts written in either Tcl or Python.

## Supported Versions
Inductiva stays up to date with the latest versions of OpenSees. Below is a list of the supported versions:

- **3.7.1** - Supports Python and Tcl scripting
- **2.5.0** - Supports Tcl scripting only

If you need to use a version not listed here, please feel free to [Contact Us](mailto:support@inductiva.ai).
We’ll be happy to accommodate your request!

## Container Images
Each version of OpenSees in the Inductiva API has its own publicly available container image, 
so you can also use it to run simulations. These images are hosted in our Docker Hub repository, 
[Kutu](https://hub.docker.com/r/inductiva/kutu/tags?name=opensees), and you can find the 
Dockerfile details for each version [here](https://github.com/inductiva/kutu/tree/main/simulators/opensees).