# Generalizing the Use Case
Inductiva enables the parallel execution of numerous simulations, making it ideal for conducting parameter studies at scale. For example, suppose you are interested in analyzing how varying a specific input parameter influences the simulation results. you might want to analyze how varying a specific input parameter affects the simulation results.

As a demonstration, consider studying the impact of changes in the reference free-stream velocity, 
`flowVelocity`, defined in the `initialConditions` file of the OpenFOAM case. This parameter represents the magnitude of the incoming flow velocity used for nondimensionalizing force coefficients. In the current configuration, `flowVelocity` is set to 20 m/s on the x axis:

```
/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Version:  8
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

flowVelocity         (20 0 0);
pressure             0;
turbulentKE          0.24;
turbulentOmega       1.78;

// ************************************************************************* //

```

We will use Inductiva’s templating system to generalize the base simulation by randomly sampling values 
for `flowVelocity` (the wind speed). This allows efficient generation of a broad set of simulation conditions.

## Parametrize the system/forceCoeffs file
Inductiva allows you to convert fixed parameters in simulation configuration files into variables that can be programmatically controlled via Python scripting.

Instead of hardcoding the wind speed in `openfoam-input-example/0/include/initialConditions` 
like this:

```
...
    flowVelocity         (20 0 0);
...
```

you can update the file to include a variable placeholder:

```
...
    flowVelocity         ({{ wind_speed }} 0 0);
...
```

After this change, save the file as `initialConditions.jinja` (note the `.jinja` extension). This tells Inductiva’s templating engine to replace the `wind_speed` variable with the appropriate value from Python.

## Code Overview
The script below shows how we can now set the value of the `flowVelocity` parameter from Python and run a variation of the original simulation:

```python
import inductiva

# Allocate cloud machine
cloud_machine = inductiva.resources.MachineGroup(
   provider="GCP",
   machine_type="c2d-highcpu-32",
   spot=True
)

wind_speed = 15

print(f"Preparing files for OpenFOAM with "\
      f"wind_speed = {wind_speed}")
target_dir = f"variations/wind_speed_{wind_speed}"

# This is where we make the substitution and set
# the value of the parameters in the modified config file
inductiva.TemplateManager.render_dir(
   source_dir="openfoam-input-example/",
   target_dir=target_dir,
   overwrite=True,
   wind_speed=wind_speed)

# Initialize OpenFOAM simulator
openfoam = inductiva.simulators.OpenFOAM()

# Run simulation
task = openfoam.run(
   input_dir=target_dir,
   shell_script="./Allrun",
   on=cloud_machine,
   resubmit_on_preemption=True,
)

task.wait()
cloud_machine.terminate()
```

That's it!

First, we set the `wind_speed` values. Then, the `render_dir` method generates the input file by replacing 
the placeholder variables in the `forceCoeffs.jinja` template with the corresponding values. The generated 
files are saved in the `target_dir` (`variations/wind_speed_{wind_speed}`).

> For more details on managing template files, check out the `TemplateManager` [documentation](https://docs.inductiva.ai/en/latest/intro_to_api/templating.html).

Since Inductiva’s API is Python-based, you can easily loop over different sampled values of `wind_speed` to run multiple simulations in parallel and generate a dataset. That’s exactly what we’ll do in the next section.