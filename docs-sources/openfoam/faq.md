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

Here's an example of how to set commands manually:

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

<br>

## 3. How do I run OpenFOAM commands in parallel?
Running OpenFOAM commands in parallel involves configuring your parallel execution settings and then running 
commands either **via a shell script** or **individually**.

### Using a shell script
You can run OpenFOAM commands in parallel within a script by either:

* Using the helper command: `runParallel <command>`
* Running with MPI directly: `mpirun -np 4 <command>`

### Running individual commands
If you're running commands manually (not via a script), you can still execute them in parallel by either:

* Using the helper: `runParallel <command>`
* Adding the `-parallel` flag directly to the command, for example: `simpleFoam -parallel`

When a command includes the `-parallel` flag, it is automatically recognized as parallel execution. 

The following two commands are equivalent:

```bash
runParallel simpleFoam
simpleFoam -parallel
```

Internally, both run as:

```bash
mpirun -np <num_processes> simpleFoam -parallel
```

### Configuring parallel execution
Before running your simulation, set the parallel settings in your machine configuration. For example:

```python
cloud_machine = inductiva.resources.MachineGroup(
    provider="GCP",
    machine_type="c2d-standard-112",
    np=5,                  # Number of processes (default: max threads)
    use_hwthread_cpus=False, # Use hyperthreading (default: True)
    mpi_version="4.1.6",     # MPI version (default: 4.1.6)
    spot=True
)
```

With this configuration, a command such as `simpleFoam -parallel` will be executed as:

```bash
mpirun -np 5 --use-hwthread-cpus simpleFoam -parallel
```

This runs the simulation using 5 processes, disables hyperthreading, and uses MPI version 4.1.6.

> **Note**: These parallel settings (number of processes, hyperthreading, and MPI version) apply only when using the `-parallel` flag. 
If you use `runParallel`, OpenFOAM manages parallelism internally and ignores these settings.

<br>


## 4. Do your simulations work on your laptop but not on Inductiva?

If your simulation runs successfully on your local machine but fails on Inductiva,
it’s likely due to a version mismatch.

Make sure you're using the same software version and distribution
(e.g., ESI or Foundation) on both your machine and Inductiva. If the version
you need isn’t available on Inductiva yet,
[contact us](mailto:support@inductiva.ai) — we’ll be happy to add it as soon as possible.

<br>

## 5. My simulation fails, but the console shows Success. What’s going on?

This usually means your `Allrun` script is not properly reporting errors.

Inductiva checks whether a simulation succeeded by looking at the script’s
**exit code**. A non-zero code means failure. However, in your case, a command
like `runApplication blockMesh` might be failing, but the script continues to
execute, often ending with a command like `echo "Simulation Complete!"`,
which returns a success code (`0`).

Because the final command succeeds, the whole script appears successful, even
though a key step failed.

### How to fix it

Add this line to the top of your `Allrun` script:

```bash
set -e
```

This will make the script stop immediately if any command fails, and it will
correctly return a non-zero exit code. That way, Inductiva can detect the
failure and show it as such.

<br>
<br>

Still can't find what you're looking for? [Contact Us](mailto:support@inductiva.ai)
