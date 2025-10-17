# ⚙️ Versions and Containers

## About GROMACS
[GROMACS](https://www.gromacs.org/index.html) is an open-source software suite widely used for simulating biochemical molecules such as proteins, lipids, and nucleic acids, but its speed and efficiency make it suitable for a range of non-biological systems, including polymers and fluid dynamics.

## Supported Versions
Inductiva stays up to date with the latest versions of GROMACS. Below is a list of the supported versions, along with their respective release dates:

- **2025.1** (Mar., 2025)
- **2025.0** (Feb., 2025)
- **2022.2** (Jun., 2022) 

> 📌 All versions above include GPU support.

To use a specific GROMACS version, simply specify it when initializing the simulator:

```python
gromacs = inductiva.simulators.GROMACS( \
    version="2025.1")
```

If you need to use a version not listed here, please feel free to [Contact Us](mailto:support@inductiva.ai).
We’ll be happy to accommodate your request!

## Container Images
Each version of GROMACS in the Inductiva API has its own publicly available container image, 
so you can also use it to run simulations. These images are hosted in our Docker Hub repository, 
[Kutu](https://hub.docker.com/r/inductiva/kutu/tags?name=gromacs), and you can find the 
Dockerfile details for each version [here](https://github.com/inductiva/kutu/tree/main/simulators/gromacs).