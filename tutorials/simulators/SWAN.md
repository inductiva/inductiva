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

## Available Commands

For this simulator, we provide the `command` argument, which specifies the
executable to run. The available commands are:

- **`swanrun`**: This is the default command when no specific command is
 provided. It typically generates debug files for your simulation. However,
 `swanrun` does not support MPI clusters, limiting its use to a single machine.

- **`swan.exe`**: This command is compatible with MPI clusters, allowing
 simulations to run across multiple machines. However, it may sometimes fail
 when generating debug files.

**Recommendation**: For simulations on a single machine, we recommend using
`swanrun`. If you need to run simulations on an MPI cluster, use `swan.exe`.

> **Note**: `swanrun` expects the simulation file in the `file.swn` format, while
`swan.exe` expects the simulation file with the name `INPUT`.

## Example Code
Below is an example of running a SWAN simulation via the Inductiva API:

```{literalinclude} ../../inductiva/tests/test_simulators/swan/swan.py
:language: python
```

## SWAN limitation with reused files across multiple simulations

With respect to [this](https://tutorials.inductiva.ai/how_to/reuse-files.html)
feature.

SWAN has constraints when reusing simulation outputs as inputs for subsequent
simulations. Currently, SWAN supports using only a single simulation output as
input for a new simulation.

The issue stems from the inability to merge outputs from multiple simulations
into a single folder. SWAN requires all simulation files to reside in the
execution folder. For instance, if you have simulation outputs in folders A and
B, running SWAN in folder A will only use the files in A, ignoring those in B.

For more detailed information on SWAN configuration, refer to the [official documentation](https://swanmodel.sourceforge.io/).

## Available Benchmarks for SWAN

The following benchmark is available to test **SWAN**:

* [Ring](https://benchmarks.inductiva.ai/SWAN/ring/): This example simulates 
wave propagation and interaction in a circular domain, based on the 
official SWAN documentation.
