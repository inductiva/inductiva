In this guide, we will explore how to set up and run FDS (Fire Dynamics Simulator), 
one of the built-in simulators available via the Inductiva API.

We will cover:

- Configuring the FDS input file for running simulations.
- Overview of the key models: hydrodynamics, combustion, and radiation transport.
- Example code to help you get started with simulations.

# FDS (Fire Dynamics Simulator)

Fire Dynamics Simulator (FDS), is a computational fluid dynamics (CFD) model of 
fire-driven fluid flow. FDS solves numerically a form of the Navier-Stokes
equations appropriate for low-speed (Ma1 < 0.3), thermally-driven flow with an
emphasis on smoke and heat transport from fires. The applications of FDS include
the design of smoke handling systems and sprinkler/detector activation studies,
as well as reconstruction of residential and industrial fires.

## Key Models in FDS

FDS contains a hydrodynamics model, a combustion model and a 
radiation transport model. All are configured by a single input file
`input_file.fds`.

**Hydrodynamic Model:** Solves the Navier-Stokes equations for low-speed, 
thermally-driven flows with an emphasis on smoke and heat transport.

**Combustion Model:** Uses a mixing-controlled chemical reaction for 
simulating combustion, with three lumped species: air, fuel, and products. 
Multiple reactions can also be modeled, depending on the simulation 
requirements.

Since FDS requires separate mesh setups assigned to individual processors, 
you will need to specify the number of cores for your simulation. Unlike 
other simulators in the Inductiva API, FDS does not automatically assign 
cores, so it’s important to configure this step manually.

## Example Code

```{literalinclude} ../../examples/fds/fds.py
:language: python
```
