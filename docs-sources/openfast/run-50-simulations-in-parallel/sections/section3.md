# Generalizing the Use Case
While a single simulation may run at similar speed, Inductiva allows you to scale many simulations in parallel. So let's say you want to study 
the effect of changing a certain input parameter of your simulation. For the sake of demonstration, let us assume that we want to study what happens 
when the offshore turbine in this example is installed in locations with different water depths.

In the 'Environmental Conditions' section of the `5MW_OC4Semi_WSt_WavesWN.fst` parameter file it can be seen that the WtrDpth parameter has been set to 200 m:

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

We will use the Inductiva API to run variations of this base simulation where we set the `WtrDpth` from 100 to 200 meters in 2 meter steps, so we will run 
50 simulations. 

## Parametrize the input file 5MW_OC4Semi_WSt_WavesWN.fst
Inductiva allows you to convert fixed parameters in your simulation configuration files into variables that you can set programmatically using Python scripting. 
This means that we will be able to change the `WtrDpth` defined in the `5MW_OC4Semi_WSt_WavesWN.fst` input file from a Python script before starting the simulation.

To do this you will need to edit your `5MW_OC4Semi_WSt_WavesWN.fst` from this:

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


After this small edit, you will need to save your input file as `5MW_OC4Semi_WSt_WavesWN.fst.jinja` (note the ".jinja" extension). 
This will tell Inductiva's templating engine that a Python variable called "water_depth" should be used to set the correct scalar value in `5MW_OC4Semi_WSt_WavesWN.fst`.

## Code Overview
The script below shows how we can now set the value of the `WtrDpth` parameter 
from Python and run a variation of the original simulation for a water depth of 190 meters instead of 200 meters:

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

Now the good thing about the Inductiva API is that it is a Python API, so you can literally just do a for loop to iterate over the whole 
range of values for `WtrDpth`. That's what we'll do next.