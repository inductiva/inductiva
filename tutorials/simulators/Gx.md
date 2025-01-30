In this guide, we will walk you through Gx setup, one of the 
built-in simulators available via the Inductiva API.

We will cover:

- Setting up Gx for use with our API.
- Example code to help you get started with simulations.

# Gx

[Gx](https://bitbucket.org/gyrokinetics/gx/src/gx/) is a high-performance code
for solving the nonlinear gyrokinetic system governing low-frequency turbulence
in magnetized plasmas. It employs Fourier-Hermite-Laguerre spectral methods,
with a key advantage being its Hermite-Laguerre velocity discretization. This
approach enables GX to seamlessly transition between coarse gyrofluid-like
resolutions and finer, more detailed gyrokinetic resolutions, providing
flexibility and accuracy in plasma turbulence simulations.


We currently have the following Gx version available:
- v11-2024 : This represents the source code from the Gx repository as of November 2024.

> **Note**: The Gx simulator only runs on VM's with a GPU available (e.g., `g2-standard-4`).

## Example Code

```{literalinclude} ../../examples/gx/gx.py
:language: python
```