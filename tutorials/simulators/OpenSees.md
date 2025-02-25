In this guide, we will walk you through setting up and running OpenSees simulations
using the Inductiva API.

We will cover an example code to help you get started with simulations.

# OpenSees

[OpenSees (Open System for Earthquake Engineering Simulation)](https://opensees.berkeley.edu/)
is an open-source software framework for simulating the response of structural
and geotechnical systems subjected to earthquakes and other dynamic loads.

Developed by the Pacific Earthquake Engineering Research (PEER) Center, OpenSees
provides a flexible and extensible platform with scripting capabilities through
Tcl and Python. It supports nonlinear analysis, finite element modeling, and
advanced material behavior, making it widely used in academic research and
engineering practice for performance-based seismic design and analysis.

## Supported Versions  
We currently support the following OpenSees versions:  
- **v2.5.0** – Supports Tcl scripting only.  
- **v3.7.1** – Supports both Python and Tcl.

## Example Code

In the following example, we demonstrate how to run an OpenFAST simulation 
using Inductiva’s cloud infrastructure. 

```{literalinclude} ../../inductiva/tests/test_simulators/opensees/opensees.py
:language: python
```

The following example demonstrates how to run an OpenSees simulation using MPI and Python.  

To use Tcl instead, simply set the `interface` argument to `tcl`.  

If you want to run a Python simulation without MPI, remove the `n_vcpus` argument
and explicitly set `use_hwthread=False`. Both of these are MPI-related arguments,
and specifying any MPI argument will result in the simulation running with MPI.