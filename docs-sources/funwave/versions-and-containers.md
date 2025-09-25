# ⚙️ Versions and Containers

## About FUNWAVE
[FUNWAVE](https://fengyanshi.github.io/build/html/index.html) is a numerical model
designed to simulate the propagation and transformation of surface gravity waves
in coastal and nearshore environments. Based on the fully nonlinear Boussinesq equations, 
it employs a Total Variation Diminishing (TVD) shock-capturing scheme, enabling it to accurately 
model highly nonlinear processes such as wave breaking, run-up, and inundation.

Thanks to its ability to capture both dispersive and nonlinear effects over complex bathymetry, 
FUNWAVE is widely used in the study of tsunamis, storm surges, and coastal wave dynamics.

## Supported Versions
Inductiva stays up to date with the latest versions of FUNWAVE. Below is a list of the supported versions, along with their respective release dates:

- **3.6** (May, 2023) 

If you need to use a version not listed here, please feel free to [Contact Us](mailto:support@inductiva.ai).
We’ll be happy to accommodate your request!

## Container Images
Each version of FUNWAVE in the Inductiva API has its own publicly available container image, 
so you can also use it to run simulations. These images are hosted in our Docker Hub repository, 
[Kutu](https://hub.docker.com/r/inductiva/kutu/tags?name=funwave), and you can find the 
Dockerfile details for each version [here](https://github.com/inductiva/kutu/tree/main/simulators/funwave).