In this guide, we will walk you through setting up and running GX simulations
using the Inductiva API.

We will cover:

- Configuring GX simulations with the appropriate input directories.
- Example code to help you get started with your simulations.

# GX

[GX](https://bitbucket.org/gyrokinetics/gx/src/gx/) is a GPU-native model for
solving the nonlinear gyrokinetic system for low-frequency turbulence in
magnetized plasmas, in particular tokamaks and stellarators, using
Fourier-Hermite-Laguerre spectral methods.

This software has proven ideal for fusion reactor design and optimization, as
well as for general physics research.

We currently have the following GX version available:
- v11-2024 : This represents the source code from the GX repository as of November 2024.

> **Note**: The GX simulator only runs on VM's with a GPU available (e.g., `g2-standard-4`).

## Example Code

In this example, we set up a linear ion-temperature-gradient (ITG) instability
calculation using W7-X stellarator geometry and adiabatic electrons. 
This use case is presented in the official GX [documentation](https://gx.readthedocs.io/en/latest/index.html),
which you can visit for more details on the simulator's features and configurations.

```{literalinclude} ../../examples/gx/gx.py
:language: python
```

The current example is divided into four steps:
   - Configure the Machine Type: In this step, we define the machine type.
   - Download Input Files: In this step, we retrieve the input files from the 
   Inductiva bucket.
   - Pick the Simulator: We select the simulator we want to use. In this case, GX.
   - Run the Simulation: In the final step, we run the simulation using the run
   method, specifying simulation configuration file.

The last three lines handle post-simulation tasks: waiting for the simulation to
finish, terminating the machine, and downloading the outputs, in that order.