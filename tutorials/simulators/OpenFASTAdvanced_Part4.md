---
orphan: true
---

# Running 40 Simulations in parallel - Templating

## Recap: Running a Single OpenFAST Simulation  

In the last part of this tutorial, we explored how to set up and run a single
OpenFAST simulation using Inductiva. We covered the necessary file preparations,
built the required DLL, and executed the simulation on a cost-effective cloud
machine. While running a single case on the cloud may not always be the best
option due to CPU clock speed limitations, the real advantage of Inductiva lies
in its ability to scale simulations effortlessly. Now, we'll leverage
this power to run hundreds of OpenFAST simulations in parallel, drastically
reducing total computation time.

## From 1 to 40
Inductiva does not help run one OpenFast simulation faster, but
it helps you run many simulations in parallel. So, let's assume
that you need to study the impact of changing a certain 
input paramter of your simulation. For the sake of demonstration,
let us assume that we want to study what happens when the off-shore
turbine of this example is installed in locations with different 
water depths. 

In the "Enviromental Conditions" section of the `5MW_OC4Semi_WSt_WavesWN.fst`
parameter file one can see that parameter WtrDpth has been set to 200 m:


```
---------------------- ENVIRONMENTAL CONDITIONS --------------------------------
    9.80665   Gravity         - Gravitational acceleration (m/s^2)
      1.225   AirDens         - Air density (kg/m^3)
       1025   WtrDens         - Water density (kg/m^3)
  1.464E-05   KinVisc         - Kinematic viscosity of working fluid (m^2/s)
        335   SpdSound        - Speed of sound in working fluid (m/s)
     103500   Patm            - Atmospheric pressure (Pa) [used only for an MHK turbine cavitation check]
       1700   Pvap            - Vapour pressure of working fluid (Pa) [used only for an MHK turbine cavitation check]
        200   WtrDpth         - Water depth (m)
          0   MSL2SWL         - Offset between still-water level and mean sea level (m) [positive upward]
```

We are going to use Inductiva API to run variations of this
base simulation where we set the WtrDpth from 180 to 220 meter,
at 1 meter steps. This means we are going to run 40 simulations.
But we are going to run these 40 simulation in parallel.

### Parametrize the input file `5MW_OC4Semi_WSt_WavesWN.fst`

Inductiva lets you transform fixed parameters in your simulation
configuration files into variables that you can set programmatically
via Python scripting. That is, we will be able to change the `WtrDpth`
defined in the `5MW_OC4Semi_WSt_WavesWN.fst` input file from a Python script
before starting the simulation.

To do that you need to edit your `5MW_OC4Semi_WSt_WavesWN.fst` from this:

```
---------------------- ENVIRONMENTAL CONDITIONS --------------------------------
...
        200   WtrDpth         - Water depth (m)
...
```

To this:

```
---------------------- ENVIRONMENTAL CONDITIONS --------------------------------
...
        {{ water_depth }}   WtrDpth         - Water depth (m)
...
```

After doing this small edit you will have to save your input file with the
following name `5MW_OC4Semi_WSt_WavesWN.fst.jinja`
(watch the ".jinja" extension). This will let Inductiva's templating 
engine know that a Python variable named "water_depth" should be 
used to set the right scalar value in `5MW_OC4Semi_WSt_WavesWN.fst`.

How is this done in practice? It's very easy. The script below
shows how we can now set the value of the WtrDpth parameter from
Python, and run a variation of the original simulation for a 
water depth of 190 meters, instaed of 200 meters:

```python
import inductiva

# Allocate cloud machine
cloud_machine = inductiva.resources.MachineGroup(
    provider="GCP",
    machine_type="n2-highcpu-2",
    spot=True
)

water_depth = 190

print(f"Preparing files for depth = {water_depth}")
target_dir = f"variations/params_for_depth_{water_depth}"

# This is where we make the substitution and set
# the value of WtrDpth in the modified config file
# 5MW_OC4Semi_WSt_WavesWN.fst.jinja
inductiva.TemplateManager.render_dir(
    source_dir="openfast-5MW_OC4Semi_WSt_WavesWN",
    target_dir=target_dir,
    overwrite=True,
    water_depth=water_depth)

# Initialize OpenFAST simulator
openfast = inductiva.simulators.OpenFAST(
    version="4.0.2")

task = openfast.run(
    input_dir=target_dir,
    commands=[
        "openfast 5MW_OC4Semi_WSt_WavesWN/"
        "5MW_OC4Semi_WSt_WavesWN.fst"],
    on=cloud_machine)

task.wait()

cloud_machine.terminate()
```

That's it!

Now, the good thing about Inductiva API is that it is a
Python API, so you can literally just do a for loop to
iterate over all the range of values for WtrDpth. That's 
what we are going to do next.

[Running 40 Simulations](OpenFASTAdvanced_Part5.md)