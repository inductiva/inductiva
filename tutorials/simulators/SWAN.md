In this guide, we will walk you through setting up and running SWAN simulations 
using the Inductiva API. 

We will cover:

- Configuring SWAN simulations using the Inductiva API.
- Example code to help you get started with simulations.
- Available benchmarks to test SWAN’s capabilities.

# SWAN (Simulating WAves Nearshore)

[SWAN](https://swanmodel.sourceforge.io/) is a simulator used to obtain 
realistic estimates of **wave parameters** in **coastal areas**, **lakes**, and 
**estuaries** based on wind, sea floor, and current conditions. SWAN is 
widely applied in coastal engineering, environmental assessments, and 
maritime safety studies.

The SWAN simulator is typically configured using a single `.swn` file, 
along with additional files that define the **domain**, **sea floor**, and **input conditions**. These files should be organized in an input directory, which 
will be passed to the simulator.

## Example Code
Below is an example of running a SWAN simulation via the Inductiva API:

```{literalinclude} ../../examples/swan/swan.py
:language: python
```

For more detailed information on SWAN configuration, refer to the [official documentation](https://swanmodel.sourceforge.io/).

## Available Benchmarks for SWAN

The following benchmark is available to test **SWAN**:

* [Ring](https://benchmarks.inductiva.ai/SWAN/ring/): This example simulates 
wave propagation and interaction in a circular domain, based on the 
official SWAN documentation.
