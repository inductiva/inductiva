In this guide, we will walk you through setting up and running SCHISM simulations using the Inductiva API. 

We will cover:

- Configuring SCHISM simulations using the Inductiva API.
- Example code to help you get started with simulations.
- Available benchmarks to test SCHISM’s capabilities in both 
2D and 3D environments.

# SCHISM

[SCHISM](https://ccrm.vims.edu/schismweb/) is an open-source 3D baroclinic 
circulation model designed to simulate hydrodynamics across creek-lake-river-estuary-shelf-ocean scales. It excels at modeling 
complex water circulation systems with varying salinity and temperature, 
making it a powerful tool for coastal, estuarine, and oceanographic 
research.

## Example Code

Below is an example of running a SCHISM simulation via the Inductiva API:

```{literalinclude} ../../examples/schism/schism.py
:language: python
```

## Available Benchmarks for SCHISM

The following benchmarks are available to test SCHISM’s performance:

* [Test_Inun_NWaves_2D](https://benchmarks.inductiva.ai/SCHISM/schism/):
The `Test_Inun_NWaves_2D` example from the SCHISM official test suite.
* [Test_Inun_NWaves_3D](https://benchmarks.inductiva.ai/SCHISM/schism_3d/):
The `Test_Inun_NWaves_3D` example from the SCHISM official test suite.
