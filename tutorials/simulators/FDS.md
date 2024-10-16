# FDS

Fire Dynamics Simulator (FDS), is a computational fluid dynamics (CFD) model of 
fire-driven fluid flow. FDS solves numerically a form of the Navier-Stokes
equations appropriate for low-speed (Ma1 < 0.3), thermally-driven flow with an
emphasis on smoke and heat transport from fires. The applications of FDS include
the design of smoke handling systems and sprinkler/detector activation studies,
as well as reconstruction of residential and industrial fires.

FDS conteplates in its model an hydrodynamics model, a combustion model and a 
radiation transport model. All are configured by a single input file
`input_file.fds`.

**Hydrodynamic Model:** solves numerically a form of the Navier-Stokes equations 
appropriate for lowspeed, thermally-driven flow with an emphasis on smoke and
heat  transport from fires. 

**Combustion Model:** For most applications, FDS uses a single step,
mixing-controlled chemical reaction which uses three lumped species (a species
representing a group of species). These lumped species are air, fuel, and
products. By default the last two lumped species are explicitly computed.
Options are available to include multiple reactions and reactions that are not
necessarily mixing-controlled.

Since FDS requires the setting of separate meshes that are attributed to the 
processors, the user needs to provide the number of cores required to run the 
simulation. This steps out of the usual structure of the simulators in Inductiva 
API, where the number of cores is automatically set.

## Example

```{literalinclude} ../../examples/fds/fds.py
:language: python
```
