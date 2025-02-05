In this guide, we will walk you through setting up and running SWASH using 
the Inductiva API. 

We will cover:

- Configuring SWASH simulations using the Inductiva API, and all 
supported versions of the simulator.
- Example code to help you get started with simulations.
- Available benchmarks to test SWASH’s capabilities.

# SWASH (Simulating WAves and Surf Hydrodynamics)

SWASH is a simulator designed to solve **shallow water equations**, used for 
simulating **waves**, **currents**, and **tidal flows** in coastal waters, harbors, 
and around coastal structures. SWASH is ideal for modeling long waves, 
tidal inlets, and rapidly varied flows in nearshore regions.

SWASH simulations are configured using a `.sws` file, with additional files 
such as **bathymetry** (with a `.bot` extension) to define the ocean floor and 
domain. These files should be placed in an input directory for the 
simulation.

## Example Code

Below is an example of running a SWASH simulation via the Inductiva API:

```{literalinclude} ../../inductiva/tests/test_simulators/swash/swash.py
:language: python
```

## Supported Versions

We currently support the following versions of SWASH:

- **9.01A** (Apr, 2023)
- **10.01A** (Apr, 2024)
- **10.05** (May, 2024)

To list all available versions of SWASH (or other simulators), you can 
use the `inductiva simulators list` CLI command.

## Available Benchmarks for SWASH

The following benchmark is available for SWASH:

* [Three-Dimensional Currents](https://benchmarks.inductiva.ai/SWASH/SWASH_Currents/):
This benchmark replicates the S1 simulation from the paper *“Modeled Three-Dimensional Currents and Eddies on an Alongshore-Variable Barred Beach”*, authored by Christine M. Baker et al. (2021).

## What to Read Next


You may also want to check the following blog post, where we illustrate a 
practical use of SWASH:

 * [Scaling coastal engineering projects with Inductiva API](https://inductiva.ai/blog/article/scaling-coastal-engineering-projects-inductiva-api)
