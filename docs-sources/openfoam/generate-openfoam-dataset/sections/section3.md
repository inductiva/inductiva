# Generalizing the Use Case
Inductiva enables the parallel execution of numerous simulations, making it ideal for conducting parameter studies at scale. For example, suppose you are interested in analyzing how varying a specific input parameter influences the simulation results.


As a demonstration, consider studying the impact of changes in the reference free-stream velocity, denoted as `flowVelocity`, which appears in the `initialConditions` file of the OpenFOAM case. This parameter represents the magnitude of the incoming flow velocity used for nondimensionalizing force coefficients. In the current configuration, flowVelocity is set to 20 m/s on the x axis:

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

We will use Inductiva’s templating system to generalize the base simulation by randomly sampling values for `flowVelocity`, representing the wind speed in the OpenFOAM case. This allows generating a broad set of simulation conditions efficiently.

## Parametrize the system/forceCoeffs file
Inductiva allows you to convert fixed parameters in your simulation configuration files into variables that can be programmatically controlled via Python scripting.

This means that instead of hardcoding the wind speed in `openfoam-input-example/0/include/initialConditions` as:
```
...
    flowVelocity         (20 0 0);
...
```


you can update the file to include a variable placeholder such as:

```
...
    flowVelocity         ({{ wind_speed }} 0 0);
...
```

After this small edit, you will need to save your input file as `initialConditions.jinja` (note the ".jinja" extension). 
This will tell Inductiva's templating engine that Python variable called "wind_speed" should be used to set the correct x axis flow velocity.

## Code Overview
The script below shows how we can now set the value of the `flowVelocity` parameter from Python:

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

First, we set the values for `wind_speed`. Then, we generate the input file using the `render_dir` method, which replaces the placeholder variables in our `forceCoeffs.jinja` template with the corresponding values. 
The generated input files are then stored in the folder, `target_dir` (`variations/wind_speed_{wind_speed}`).
 For more details on managing template files, check out the `TemplateManager` [documentation](https://docs.inductiva.ai/en/latest/intro_to_api/templating.html).

The good thing about the Inductiva API is that it’s a Python API, so you can simply use a for loop with random sampling for `wind_speed` to run multiple simulations in parallel and generate a dataset. That’s exactly what we’ll do in the next section.
