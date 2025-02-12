In this guide, we will walk you through setting up and running OpenFAST, 
an open-source engineering toolset developed by the National Renewable 
Energy Laboratory (NREL) for simulating wind turbine dynamics.

We will cover:

- Configuring OpenFAST and FAST.Farm for wind turbine and wind farm simulations.
- Allowed Commands for running OpenFAST.
- Example code to help you get started with simulations.
- Available benchmarks to help you evaluate the performance of OpenFAST.

# OpenFAST

[OpenFAST](https://www.nrel.gov/wind/nwtc/openfast.html) provides a 
modular framework for designing and analyzing wind turbines, allowing for 
detailed modeling of components such as blades, towers, and control systems.

OpenFAST comes integrated with FAST.Farm, an extension designed for 
large-scale wind farm simulations. FAST.Farm models complex interactions 
between turbines, accounting for factors like wake effects, turbulence, 
and terrain. FAST.Farm's capabilities are indispensable for the design 
and planning of wind energy projects, unlocking the full potential of 
wind resources for renewable energy generation.

Both OpenFAST and FAST.Farm are compiled with OpenMP to enable parallel processing, enhancing the performance and scalability of your simulations.

## Allowed Commands

The following commands are available for running OpenFAST and its modules:

- `aerodyn_driver`: Performs **aerodynamic analysis** for airflow dynamics 
around wind turbine blades.
- `beamdyn_driver`: Focuses on **structural analysis** for the dynamic 
response of wind turbine blades and towers.
- `feam_driver`: Executes **finite element analysis (FEA)**, focusing on 
detailed structural modeling and simulation.
- `hydrodyn_driver`: Simulates **hydrodynamic interactions** between turbine 
support structures and water bodies.
- `inflowwind_driver`: Simulates **atmospheric conditions**, such as wind 
inflow, and their effects on wind turbine performance.
- `moordyn_driver`: Focuses on **mooring dynamics** for floating wind turbines, 
analyzing the behavior of mooring systems.
- `openfast`: The main **OpenFAST** executable for integrating and running 
comprehensive wind turbine simulations.
- `orca_driver`: Enables **aero-elastic and hydrodynamic coupled simulations**,
useful for analyzing turbine behavior in complex marine environments.
- `servodyn_driver`: Simulates **control system dynamics**, modeling how 
wind turbine control systems respond to changing conditions.
- `subdyn_driver`: Focuses on **substructure dynamics**, analyzing the interaction 
between turbine components and support structures.
- `turbsim`: A standalone tool that generates **atmospheric turbulence** data 
for wind turbine simulations.
- `unsteadyaero_driver`: Used for **unsteady aerodynamics analysis**, focusing 
on time-varying airflow around turbine blades.
- `FAST.Farm`: Specialized for **wind farm simulations**, accounting for 
turbine interactions, wake effects, and turbulence.

## Example Code

In the following example, we demonstrate how to run an OpenFAST simulation 
using Inductiva’s cloud infrastructure. 

```{literalinclude} ../../inductiva/tests/test_simulators/openfast/openfast.py
:language: python
```

## Inductiva OpenFAST Benchmarks

The following benchmarks are available to help you evaluate the performance of 
the OpenFAST suite on Inductiva’s infrastructure:

* [5MW Land](https://benchmarks.inductiva.ai/OpenFAST/OpenFAST_Land/):
A simulation of a land-based NREL 5-MW turbine using **BeamDyn** as the 
structural module. It simulates 20 seconds with a time step size of 0.001.
* [5MW OC4](https://benchmarks.inductiva.ai/OpenFAST/OpenFAST_OC4/): Simulates 
an offshore, fixed-bottom NREL 5-MW turbine, with a focus on **HydroDyn** wave dynamics.
* [FAST.Farm](https://benchmarks.inductiva.ai/OpenFAST/OpenFAST_FAST.Farm/): 
A tool for forecasting the performance and loads of wind turbines within 
a wind farm, running simulations in parallel with **OpenFAST**.

