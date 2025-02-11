In this guide, we will walk you through AMR-Wind’s setup, one of the 
built-in simulators available via the Inductiva API.

We will cover:

- Configuring Amr-Wind simulations using the Inductiva API.
- Example code to help you get started with simulations.

# Amr-Wind

[AMR-Wind](https://github.com/Exawind/amr-wind) is a highly parallel, 
block-structured adaptive-mesh solver specifically designed for simulating 
**incompressible** flows in wind turbines and wind farms. It is built on 
**incflo** and optimized for wind-related computations.

AMR-Wind utilizes the **AMReX library**, which provides the necessary mesh 
data structures, adaptive mesh refinement (AMR) capabilities, and linear 
solvers needed to solve the governing equations efficiently.

We are currently running **AMR-Wind version 1.4.0** for **CPU**, and support for 
**GPU** is coming soon. The solver is compiled with **OpenMPI 4.1.2** to ensure 
robust performance in parallel simulations.

## Amr-Wind Example Code

Below is an example of running an Amr-Wind simulation via the Inductiva API:

```{literalinclude} ../../inductiva/tests/test_simulators/amr_wind/amr_wind.py
:language: python
```

```
