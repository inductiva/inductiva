In this guide, we will walk you through CaNS setup, one of the 
built-in simulators available via the Inductiva API.

We will cover:

- Setting up CaNS for use with our API.
- Example code to help you get started with simulations.

# CaNS

[CaNS](https://github.com/CaNS-World/CaNS) (Canonical Navier-Stokes) is a
simulator designed for **massively-parallel numerical simulations** of fluid 
flows. CaNS specializes in simulating **incompressible, Newtonian fluid** 
**flows** and leverages a highly efficient **FFT-based solver** for the 
second-order finite-difference Poisson equation on a 3D Cartesian grid.


CaNS is currently running **version 2.3.4 for CPU** and is built using **OpenMPI 4.1.2**
and **FFTW 3.3.8**. A **GPU version** is planned for a future release to 
enhance performance.

## Example Code

```{literalinclude} ../../inductiva/tests/test_simulators/cans/cans.py
:language: python
```

**Closing Notes**: Make sure to include a designated ‘data’ folder in 
your input files to store simulation outputs.