# Versions and Containers 🛠️

## Supported Versions
Inductiva stays up to date with the latest versions of Quantum ESPRESSO. Below is a list of the supported versions, along with their respective release dates:

- **7.4.1** (Feb., 2025) for CPU
- **7.3.1** (Feb., 2024) for CPU

Inductiva has compiled Quantum ESPRESSO for two parallel execution models: **MPI** and **OpenMP**. To run the MPI version, use the standard command names (e.g., `pw.x`). For the OpenMP version, append `_openmp` to the command names (e.g., `pw_openmp.x`). 

If you need to use a version not listed here, please feel free to [Contact Us](mailto:support@inductiva.ai).
We’ll be happy to accommodate your request!

## Container Images
Each version of Quantum ESPRESSO in the Inductiva API has its own publicly available container image, 
so you can also use it to run simulations. These images are hosted in our Docker Hub repository, 
[Kutu](https://hub.docker.com/r/inductiva/kutu/tags?name=quantum-espresso), and you can find the 
Dockerfile details for each version [here](https://github.com/inductiva/kutu/tree/main/simulators/quantum-espresso).