**Find answers to commonly asked questions about OpenFOAM.**

<br>

# FAQ

## 1. What file structure is assumed for OpenFOAM simulations?
We assume the canonical OpenFOAM directory structure, which includes the following key folders:

- `time`: Contains files for particular fields, such as
initial values and boundary conditions that you must specify. For example, 
initial conditions at  t=0  are stored in the `0` directory.
- `constant`: Contains files that describe the objects in the simulation and the 
physical properties being modeled.
- `system`: Contains files that describe the simulation, including solvers, 
numerical parameters, and output files. 

<br>

## 2. Can I manually set and run OpenFOAM commands?
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

## 3. How Can I Run Parallel Commands in OpenFOAM?

OpenFOAM commands can be run in parallel in a couple of ways.

### a. If You're Using a Shell Script

You can run commands in parallel by either:

* Using the helper: `runParallel <command>`, or
* Running the command directly with MPI: `mpirun -np 4 <command>`

### b. If You're Providing a List of Commands

When submitting individual commands (not via a script), there are two main options:

* Use `runParallel <command>`, just like with scripts, or
* Use the `-parallel` flag directly with the command, for example: `simpleFoam -parallel`

If a command includes the `-parallel` flag, it will automatically be interpreted
as a parallel execution. For example, both of the following are equivalent:

```bash
runParallel simpleFoam
simpleFoam -parallel
```

Internally, this is converted to a command like:

```bash
mpirun -np <num_processes> simpleFoam -parallel
```

### Setting Up Parallelism

Parallel behavior is configured when defining the machine on which the simulation will run:

```python
cloud_machine = inductiva.resources.MachineGroup(
    provider="GCP",
    machine_type="c2d-standard-112",
    np=5, # Default is the maximum number of threads available
    use_hwthread_cpus=False, # Default is True
    mpi_version="4.1.6", # Default is 4.1.6
    spot=True
)
```

With this configuration, a command like `simpleFoam -parallel` will be run as:

```bash
mpirun -np 5 --use-hwthread-cpus simpleFoam -parallel
```

This runs the simulation using 5 processes, without hyperthreading, and using
MPI version 4.1.6.

> **Note:** Parallel settings like number of processes, MPI version, and hyperthreading only apply when using `<command> -parallel`.
> If you use `runParallel`, OpenFOAM handles the parallel execution internally and does not use the settings above.

<br>
<br>

Still can't find what you're looking for? [Contact Us](mailto:support@inductiva.ai)
