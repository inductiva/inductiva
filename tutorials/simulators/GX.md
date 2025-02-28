This guide will walk you through setting up and running GX simulations
using the Inductiva API.

We will cover an example code to help you get started with simulations.

# GX

[GX](https://bitbucket.org/gyrokinetics/gx/src/gx/) is a GPU-native model for
solving the non-linear gyrokinetic system for low-frequency turbulence in
magnetized plasmas, in particular tokamaks and stellarators, using
Fourier-Hermite-Laguerre spectral methods.

This software has proven ideal for fusion reactor design and optimization, as
well as for general physics research.

We currently have the following GX version available:
- **v11-2024** - This version corresponds to the source code from the GX repository as of November 2024.

> **Note**: The GX simulator only runs on VMs with a GPU available (e.g., `g2-standard-4`).

## Example Code

This example demonstrates how to run a GX simulation using a linear use case available 
in the official [GX documentation](https://gx.readthedocs.io/en/latest/LinearStell.html). 
Before you start, download the input files [here](https://bitbucket.org/gyrokinetics/gx/src/gx/benchmarks/linear/ITG_w7x/).

Here is the code required to run a GX simulation using the Inductiva API:

```{literalinclude} ../../inductiva/tests/test_simulators/gx/gx.py
:language: python
```

To adapt it for this or any other use case, simply replace `input_dir` with the path to your GX files and specify the `sim_config_filename` before running it in a Python script.

Once the simulation is complete, we terminate the machine, download the results and print a summary of the simulation.

It’s that simple!