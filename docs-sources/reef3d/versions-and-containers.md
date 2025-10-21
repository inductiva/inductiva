# ⚙️ Versions and Containers

## About REEF3D
[REEF3D](https://reef3d.wordpress.com/) is an open-source hydrodynamics framework specifically designed for coastal, marine, and hydraulic engineering applications. 
Built with a modular programming approach, it offers multiphysics solvers tailored to address fluid flow problems, such as sediment transport and floating body dynamics, as well as wave modeling.

The modular programming approach allows the framework to incorporate a range of different flow solvers which together represent all relevant length scales. Depending on the wave or flow conditions, the following optimized hydrodynamic modules are available:
- **REEF3D::CFD** solves the Navier-Stokes equations in three dimensions. For near-field simulations with a complex free surface pattern, it uses a two-phase flow approach with the level set method for interface capturing.
- **REEF3D::FNPF** is a three-dimensional fully nonlinear potential flow solver. It is massively parallelized and can be used to create large-scale phase-resolved sea states at all water depths.
- **REEF3D::SFLOW** is a depth-averaged model, solving the non-hydrostatic shallow water equations ideal for near-shore hydrodynamics and river flow.

## Supported Versions
Inductiva stays up to date with the latest versions of REEF3D. Below is a list of the supported versions, along with their respective release dates:

- **25.07** (Jul., 2025)
- **25.05** (Jun., 2025)
- **25.02** (Feb., 2025)
- **24.12** (Dec., 2024)
- **24.08** (Sep., 2024)
- **24.05** (Jun., 2024)
- **24.03** (Mar., 2024)
- **24.02** (Mar., 2024)

To use a specific REEF3D version, simply specify it when initializing the simulator:

```python
reef3d = inductiva.simulators.REEF3D( \
    version="25.07")
```

If you need to use a version not listed here, please feel free to [Contact Us](mailto:support@inductiva.ai).
We’ll be happy to accommodate your request!

## Container Images
Each version of REEF3D in the Inductiva API has its own publicly available container image, 
so you can also use it to run simulations. These images are hosted in our Docker Hub repository, 
[Kutu](https://hub.docker.com/r/inductiva/kutu/tags?name=reef3d), and you can find the 
Dockerfile details for each version [here](https://github.com/inductiva/kutu/tree/main/simulators/reef3d).