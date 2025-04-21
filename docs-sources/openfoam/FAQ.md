**Find answers to commonly asked questions about OpenFOAM.**

<br>

# FAQ

## What file structure is assumed for OpenFOAM simulations?
We assume the canonical OpenFOAM directory structure, which includes the following key folders:

- `time`: Contains files for particular fields, such as
initial values and boundary conditions that you must specify. For example, 
initial conditions at  t=0  are stored in the `0` directory.
- `constant`: Contains files that describe the objects in the simulation and the 
physical properties being modeled.
- `system`: Contains files that describe the simulation, including solvers, 
numerical parameters, and output files. 

<br>

## Can I manually set and run OpenFOAM commands?
Yes, absolutely! For greater flexibility, you can manually define and run each command step-by-step, giving you full control over your simulation.

Here’s an example of how to set commands manually:

```python
commands_single_machine = [
    "runApplication surfaceFeatures",
    "runApplication blockMesh",
    "runApplication decomposePar -copyZero",
    "runParallel snappyHexMesh -overwrite",
    "runParallel potentialFoam",
    "runParallel simpleFoam",
    "runApplication reconstructParMesh -constant",
    "runApplication reconstructPar -latestTime"
]

task = openfoam.run(
    input_dir=input_dir,
    commands=commands_single_machine,
    on=cloud_machine)
```

For more details on available commands and MPI configuration, check out the Custom Docker Images documentation.

<br>
<br>

Still can't find what you're looking for? [Contact Us](mailto:support@inductiva.ai)
