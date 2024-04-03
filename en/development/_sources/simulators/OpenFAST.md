# OpenFAST

[OpenFAST](https://www.nrel.gov/wind/nwtc/openfast.html) is an open-source 
engineering toolset developed by the National Renewable Energy Laboratory (NREL)
for simulating wind turbine dynamics. It provides a modular framework for designing
and analyzing wind turbines, allowing engineers to model various components such
as blades, towers, and control systems.

OpenFAST has been compiled with FAST.Farm, an extension tailored for large-scale
wind farm simulations. FAST.Farm facilitates the modeling of interactions between
multiple turbines within wind farms, accounting for factors like wake effects,
turbulence, and terrain complexity. FAST.Farm's capabilities are indispensable
for the design and planning of wind energy projects, unlocking the full potential
of wind resources for renewable energy generation.

Moreover, both OpenFAST and FAST.Farm were compiled with OpenMP enabled, further
optimizing their performance and scalability. This integration of parallel
processing techniques enhances the efficiency and speed of simulations,
empowering researchers and developers to tackle complex wind energy challenges
with greater precision and accuracy.

### Example


```python
import inductiva

# Set simulation input directory
input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "openfast-input-example.zip", unzip=True)

# List of commands to run
commands = [
    "openfast IEA-15-240-RWT-Monopile.fst"
]

# Initialize OpenFAST simulator
openfast = inductiva.simulators.OpenFAST()

# Run simulation 
task = openfast.run(input_dir=input_dir,
                    commands=commands)

task.wait()
task.download_outputs()

```

### Allowed Commands

- `aerodyn_driver`: Utilized for aerodynamic analysis, this binary focuses on 
airflow dynamics around wind turbine blades.
- `beamdyn_driver`: Primarily used for structural analysis, this binary focuses on
the dynamic response of wind turbine blades and towers.
- `feam_driver`: This binary is employed for finite element analysis, focusing on
detailed structural modeling and simulation.
- `hydrodyn_driver`: Specifically designed for hydrodynamic analysis, this binary
simulates the interaction between wind turbine support structures and water
bodies.
- `inflowwind_driver`: Used for inflow wind modeling, this binary simulates
atmospheric conditions and their effects on wind turbine performance.
- `moordyn_driver`: Focused on mooring dynamics, this binary analyzes the behavior
of floating wind turbines and their mooring systems.
- `openfast`: The main OpenFAST executable, orchestrating the integration and
execution of various modules for comprehensive wind turbine simulation.
- `orca_driver`: Employed for aero-elastic and hydrodynamic coupled simulations,
this binary enables advanced analysis of wind turbine behavior in complex 
environments.
- `servodyn_driver`: This binary specializes in control system analysis, simulating
the response of wind turbine control systems to varying conditions.
- `subdyn_driver`: Utilized for substructure dynamics analysis, this binary focuses
on the interaction between wind turbine components and their support structures.
- `turbsim`: A standalone tool for simulating atmospheric turbulence and its effects
on wind turbine performance.
- `unsteadyaero_driver`: This binary is used for unsteady aerodynamics analysis,
focusing on time-varying airflow around wind turbine blades.
- `FAST.Farm`: An extension of OpenFAST tailored for large-scale wind farm
simulations, facilitating the modeling of interactions between multiple turbines
within wind farms, accounting for factors like wake effects, turbulence, and
terrain complexity.