# Generalizing the Use Case
Inductiva enables the parallel execution of numerous simulations, making it ideal for conducting parameter studies at scale. For example, suppose you are interested in analyzing how varying a specific input parameter influences the simulation results.


As a demonstration, consider studying the impact of changes in the mean wind velocity at the reference height, denoted as `URef`. This parameter is defined in the Meteorological Boundary Conditions section of the `90m_12mps_twr.inp` input file. In the current configuration, `URef` is set to 12 m/s.

```
--------Meteorological Boundary Conditions-------------------
"IECKAI"      TurbModel       - Turbulence model
"1-ed3"       IECstandard     - Number of IEC 61400-x standard
"B"           IECturbc        - IEC turbulence characteristic
"NTM"         IEC_WindType    - IEC turbulence type 
"default"     ETMc            - IEC Extreme Turbulence Model "c" parameter [m/s]
"PL"          WindProfileType - Velocity profile type 
         90   RefHt           - Height of the reference velocity (URef) [m]
         12   URef            - Mean (total) velocity at the reference height [m/s]
"default"     ZJetMax         - Jet height [m]
"default"     PLExp           - Power law exponent [-]
"default"     Z0              - Surface roughness length [m]

```

We will use Inductiva templating system to generalize the base simulation by randomly sampling values for `URef`, `RandSeed1`, and `RandSeed2`. In the next section we will use this approach to generate a diverse dataset of simulations.

## Parametrize the input file 90m_12mps_twr.inp
Inductiva allows you to convert fixed parameters in your simulation configuration files into variables that you can set programmatically using Python scripting. 
This means that we will be able to change the `URef`, `RandSeed1`, and `RandSeed2` defined in the `90m_12mps_twr.fst` input file from a Python script before starting the simulation.

To do this you will need to edit your `90m_12mps_twr.inp` from this:

```
---------Runtime Options-----------------------------------
...
      13428  RandSeed1        - First random seed
   "RanLux"  RandSeed2        - Second random seed
...
--------Meteorological Boundary Conditions-------------------
...
         12  URef             - Mean (total) velocity at the reference height [m/s]
...
```


To this:


```
---------Runtime Options-----------------------------------
...
   {{ seed_1 }}  RandSeed1    - First random seed
   {{ seed_2 }}  RandSeed2    - Second random seed
...
--------Meteorological Boundary Conditions-------------------
...
   {{ URef }}    URef         - Mean (total) velocity at the reference height [m/s]
...
```




After this small edit, you will need to save your input file as `90m_12mps_twr.inp.jinja` (note the ".jinja" extension). 
This will tell Inductiva's templating engine that Python variables called "URef", "seed_1" and "seed_2" should be used to set the correct scalar value in `90m_12mps_twr.fst`.

## Code Overview
The script below shows how we can now set the value of the `URef`, `RandSeed1`, and `RandSeed2` parameters from Python:

```python
import inductiva

# Allocate cloud machine
cloud_machine = inductiva.resources.MachineGroup(
   provider="GCP",
   machine_type="n2-highcpu-2",
   spot=True
)

seed_1 = 3
seed_2 = 5
URef = 15

print(f"Preparing files for TurbSim with "\
      f"seed1 = {seed_1}, seed_2 = {seed_2}, "\
      f"URef = {URef}")
target_dir = f"variations/s1_{seed_1}/s2_{seed_2}/URef_{URef}"

# This is where we make the substitution and set
# the value of the parameters in the modified config file
# turbsim 90m_12mps_twr.inp.jinja
inductiva.TemplateManager.render_dir(
   source_dir="input_files",
   target_dir=target_dir,
   overwrite=True,
   seed_1=seed_1,
   seed_2=seed_2,
   URef=URef)

# Initialize OpenFAST simulator
turbsim = inductiva.simulators.OpenFAST(version="3.5.2")

task = openfast.run(
   input_dir=target_dir,
   commands=["turbsim 90m_12mps_twr.inp"],
   on=cloud_machine)

task.wait()
cloud_machine.terminate()
```

That's it!

First, we set the values for `seed_1`, `seed_2` and `URef`. Then, we generate the input file using the `render_dir` method, which replaces the placeholder variables in our `90m_12mps_twr.fst.jinja` template with the corresponding values. 
The generated input files are then stored in the folder, `target_dir` (`variations/s1_{seed_1}/s2_{seed_2}/URef_{URef}`).
 For more details on managing template files, check out the `TemplateManager` [documentation](https://docs.inductiva.ai/en/latest/intro_to_api/templating.html).

The good thing about the Inductiva API is that it’s a Python API, so you can simply use a for loop with random sampling for `URef`, `seed_1` and `seed_2` to run multiple simulations in parallel and generate a dataset. That’s exactly what we’ll do in the next section.
