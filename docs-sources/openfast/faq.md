**Find answers to commonly asked questions about OpenFAST.**

<br>

# FAQ

## 1. Can OpenFAST be run in parallel?
Yes, OpenFAST can be run in parallel, but only under specific conditions. To
utilize parallel processing, OpenFAST requires either OpenMP support via OLAF
(OpenFAST Large-Area Farm) or the FAST.Farm framework. These tools enable
parallel execution. If you are not using OLAF or FAST.Farm, OpenFAST will run
in a single-threaded mode and cannot take advantage of parallel computing capabilities.

*Source: [OpenFAST Documentation](https://openfast.readthedocs.io/en/main/source/user/fast.farm/Introduction.html#fast-farm-parallelization)*

<br>

## 2. What is the best machine for running OpenFAST?
The ideal machine for running OpenFAST depends on how you plan to run it. If
you're using OpenFAST with parallel processing (through OLAF or FAST.Farm),
a machine with a higher number of cores will provide better performance. However,
if you’re running OpenFAST in serial mode (without parallel support), a more
affordable machine with fewer cores will suffice, as OpenFAST won't utilize
additional cores in this mode.

<br>

## 3. Can I run multiple OpenFAST simulations in parallel?
Absolutely! Running multiple simulations in parallel is the best way to fully
utilize your available cloud resources. While OpenFAST can have limitations in
terms of parallelism for individual simulations, executing several simulations
concurrently allows you to maximize efficiency. For detailed guidance on how to
set this up, check out this [OpenFAST tutorial](https://inductiva.ai/guides/openfast/OpenFAST_advanced)
to see how you can do it.

<br>

## 4. How to Load the `libdiscon.so` library?

Some simulations require the `libdiscon.so` library. To simplify the process,
we've included this library in our simulation image. You can use it by
referencing the following path in your simulation files:  

```
/usr/lib/x86_64-linux-gnu/libdiscon.so
```

<br>

## 5. What commands are used to run OpenFAST and its modules?
Here is a list of the commands available to run OpenFAST v4.0.2 and its modules, in alphabetical order:

- `aerodisk_driver`: Performs **aerodynamic analysis** for wind turbine rotors, using an actuator disk approach
 to simplify the calculation of forces and moments
- `aerodyn_driver`: Performs **aerodynamic analysis** for airflow
  dynamics around wind turbine blades
- `beamdyn_driver`: Focuses on **structural analysis** for the dynamic
  response of wind turbine blades and towers
- `feam_driver`: Executes **finite element analysis (FEA)**, focusing
  on detailed structural modeling and simulation
- `FAST.Farm`: Specialized for **wind farm simulations**, accounting
  for turbine interactions, wake effects, and turbulence
- `hydrodyn_driver`: Simulates **hydrodynamic interactions** between
  turbine support structures and water bodies
- `inflowwind_driver`: Simulates **atmospheric conditions**, such as
  wind inflow, and their effects on wind turbine performance
- `moordyn_driver`: Focuses on **mooring dynamics** for floating wind
  turbines, analyzing the behavior of mooring systems
- `openfast`: The main **OpenFAST** executable for integrating and
  running comprehensive wind turbine simulations
- `orca_driver`: Enables **aero-elastic and hydrodynamic coupled
  simulations**, useful for analyzing turbine behavior in complex
  marine environments
- `seastate_driver`: Generates **wave field information** used by HydroDyn
- `sed_driver`: Simulates **simplified structural dynamics**, focusing on rigid rotor motion and enabling faster simulations
  with larger time steps
- `servodyn_driver`: Simulates **control system dynamics**, modeling
  how wind turbine control systems respond to changing conditions
- `subdyn_driver`: Focuses on **substructure dynamics**, analyzing the
  interaction between turbine components and support structures
- `turbsim`: A standalone tool that generates **atmospheric turbulence**
  data for wind turbine simulations
- `unsteadyaero_driver`: Used for **unsteady aerodynamics analysis**,
  focusing on time-varying airflow around turbine blades

<br>
<br>

Still can't find what you're looking for? [Contact Us](mailto:support@inductiva.ai)
