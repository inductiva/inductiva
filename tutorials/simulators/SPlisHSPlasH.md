In this guide, we will walk you through setting up and running 
SPlisHSPlasH, a Smoothed-Particle Hydrodynamics (SPH) simulator 
available via the Inductiva API.

We will cover:

- Configuring SPlisHSPlasH simulations using the Inductiva API.
- Example code to help you get started with simulations.
- Available benchmarks to test SPlisHSPlasH’s performance.

# SPlisHSPlasH

SPlisHSPlasH is a Smoothed-Particle Hydrodynamics (SPH) simulator that 
is widely used for simulating fluid dynamics across a variety of applications. 
Simulations are typically configured using a single `.json` file, which 
defines the geometry, fluid properties, boundary conditions, numerical 
parameters, and output settings. Additional geometry files may also be 
required in some cases.

## Example Code

Below is an example of running a simple SPlisHSPlasH simulation via the 
Inductiva API:

```{literalinclude} ../../inductiva/tests/test_simulators/splishsplash/splishsplash.py
:language: python
```

## Available Benchmarks for SPlisHSPlasH

The following benchmarks are available to test SPlisHSPlasH’s performance:

* [Fluid Cube S](https://benchmarks.inductiva.ai/SPlisHSPlasH/splish_splash/):
This benchmark replicates the example in [this tutorial](https://tutorials.inductiva.ai/generating-synthetic-data/synthetic-data-generation-1.html), using the deafult values.
* [Fluid Cube L](https://benchmarks.inductiva.ai/SPlisHSPlasH/splish_splash/):
This benchmark mirrors the benchmark [Fluid Cube S](https://benchmarks.inductiva.ai/SPlisHSPlasH/splish_splash/) in all aspects except for the particle radius, which has been decreased to 0.0045 and the fluid model, that was doubled in all axis.
* [Fluid Cube M](https://benchmarks.inductiva.ai/SPlisHSPlasH/splish_splash/): This benchmark mirrors the [Fluid Cube S](https://benchmarks.inductiva.ai/SPlisHSPlasH/splish_splash/) benchmark in all aspects except for the particle radius, which has been decreased to 0.0045.

## What to Read Next

You may want to explore the following tutorial, where we demonstrate how 
to generate synthetic data for training Physics-ML models using SPlisHSPlasH:

 * [Generating Synthetic Data for training Physics-ML models](https://tutorials.inductiva.ai/generating-synthetic-data/synthetic-data-generation-1.html)