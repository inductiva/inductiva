In this guide, we will walk you through setting up and running OpenFOAM 
simulations using the Inductiva API.

We will cover:

- Configuring OpenFOAM simulations with the appropriate input directories.
- Example codes to help you get started with simulations.
- An advanced example for running MB9 Micro-benchmark by ExaFOAM.

# OpenFOAM

OpenFOAM is a Finite Volume method for CFD simulations with a wide range of 
applications across several areas of engineering and science. It offers 
a broad set of features for everything from **complex fluid flows** (including 
chemical reactions, turbulence, and heat transfer) to **solid dynamics** and 
**electromagnetics**.

There are two main open-source distributions of OpenFOAM: one developed by the
[OpenFOAM foundation](https://openfoam.org/) and another by the
[ESI Group](https://www.openfoam.com/). The Inductiva API supports both,
and you can select your preferred distribution by setting the `distribution` parameter
when initializing the simulator. *By default, it uses the OpenFOAM Foundation version.*

We are assuming the canonical file structure for OpenFOAM simulations, which
includes the `time`, `constant`, and `system` directories.
- `time`: Contains files for particular fields, like
initial values and boundary conditions that you must specify. For example, 
initial conditions at  t=0  are stored in the `0` directory.
- `constant`: Contains files that describe the objects in the simulation and the 
physical properties being modeled.
- `system`: Contains files that describe the simulation, including solvers, 
numerical parameters, and output files. 

In order to run your simulation you cam simply run your `Allrun` script and you
are good to go.

## First Example: Motorbike Tutorial

In this example, we demonstrate how to run the [motorbike tutorial](https://github.com/OpenFOAM/OpenFOAM-8/tree/master/tutorials/incompressible/simpleFoam/motorBike) 
tutorial using the OpenFOAM Foundation distribution.

```{literalinclude} ../../examples/openfoam_foundation/openfoam_foundation.py
:language: python
```

### Example Code - ESI Distribution

To run the sample simulation above, simply download the
`openfoam-esi-input-example.zip` file, select the correct distribution by
using `inductiva.simulators.OpenFOAM(distribution="esi")`, and run the respetive
`Allrun` script.

## Advanced Tutorial: Running the MB9 Micro-benchmark from ExaFOAM

This guide walks you through running a complex OpenFOAM simulation using the
**MB9 micro-benchmark** from [ExaFOAM](https://exafoam.eu/benchmarks/). This
benchmark simulates a high-lift aircraft configuration, ideal for studying
near-wall turbulence using **wall-modeled Large Eddy Simulation (WMLES)**.

### Objective

We'll run this simulation on:
1. **Single 360 vCPU Machine**.
2. **MPI Cluster** using two 360 vCPU machines to improve performance.

### Prerequisites

1. **Download Input Files**: Get the input files from the
[repository](https://develop.openfoam.com/committees/hpc/-/tree/develop/compressible/rhoPimpleFoam/LES/highLiftConfiguration)
and place them in a folder named `highLiftConfiguration`.

   **Directory Structure**:
   ```bash
   ls -lasgo highLiftConfiguration
   total 104
   0 drwxrwxr-x@ 14     448 Sep 23 11:44 .
   0 drwx------@ 19     608 Sep 23 11:49 ..
   0 drwxrwxr-x@ 12     384 Sep 20 09:43 0.orig
   8 -rwxr-xr-x@  1     626 Jun 13 10:09 Allclean
   16 -rwxr-xr-x@  1    6998 Jun 13 10:09 Allrun
   8 -rw-rw-r--@  1     991 Jun 13 10:09 COPYING
   48 -rw-rw-r--@  1   21547 Jun 13 10:09 README.md
   0 -rw-rw-r--@  1       0 Jun 13 10:09 case.foam
   0 drwxrwxr-x@  6     192 Sep 20 09:43 constant
   0 drwxrwxr-x@ 13     416 Jun 13 10:09 figures
   0 drwxrwxr-x@ 28     896 Sep 20 09:43 system
   24 -rw-rw-r--@  1   11399 Jun 13 10:09 thumbnail.png
   ```

### Overview

Here’s the code you'll be working on as we progress through the tutorial. Don’t
worry if it doesn’t all make sense right now; everything will become clearer
in the upcoming steps.

```python
import inductiva

machine_group = inductiva.resources.MachineGroup(
            machine_type="c3d-highcpu-360",
            spot=True)
machine_group.start()

input_dir = "/path/to/highLiftConfiguration"

#Choose your simulator
openfoam = inductiva.simulators.OpenFOAM(distribution="esi")

task = openfoam.run(
            input_dir=input_dir,
            bash_script="./Allrun",
            n_vcpus=180,
            use_hwthread=True,
            on=machine_group)

task.wait()
task.download_outputs()
machine_group.terminate()

task.print_summary()

```

### Step 1: Adjust Simulation Parameters

For a faster simulation, modify the following parameters in the case definition
file (`system/include/caseDefinition`):

- **Time Step (`dt`)**: Set to 0.00002.
- **Start Time (`initTime`)**: 0.10
- **End Time (`finalTime`)**: 0.30

### Step 2: Running the Simulation

#### a. Configure and Start Machine

1. **Pick your machine**:
    ```python
    import inductiva
    machine_group = inductiva.resources.MachineGroup(machine_type="c3d-highcpu-360", spot=True)
    ```
    **Note**: `spot` machines are a lot cheaper but can be terminated by the
    provider if needed.

2. **Start your machine**
    ```python
    machine_group.start()
    ```

#### b. Simulation inputs
1. **Specify Simulation Directory**:
Let's start by defining a variable that points to the `highLiftConfiguration`
folder where all your simulation files are located.

   ```python
   input_dir = "/path/to/highLiftConfiguration"
   ```
2. **Bash script**:
To run this simulation, simply execute the `Allrun` file by pointing to it like
so:
   ```python
    bash_script = "./Allrun"
   ```

#### c. Run your simulation

1. **Run the simulation**:
We now have all we need to run our simulation.
   ```python
   #Choose your simulator
   openfoam = inductiva.simulators.OpenFOAM(distribution="esi")

   task = openfoam.run(
               input_dir=input_dir,
               bash_script="./Allrun",
               on=machine_group)
   ```

2. **Wait and Download Outputs**:
That is it. Our simulation is now running on the cloud. We can `wait` for the
simulation to be over, or we can turn our computer off go for a coffe (☕️).
   ```python
   task.wait()
   task.download_outputs()
   ```

3. **Terminate Machine**:
Once our simulation is over we can/should terminate our machine to save on costs.
If you forget, dont worry we got your back. By default, a machine will be
automaticly terminated if no simulation runs on it for 30 minutes.

   ```python
   machine_group.terminate()
   ```

4. **Check your simulation summary**:
Now that our simulation has finished we can print a summary of said simulation.
This includes information about the execution times, outputs generated and
much more.
    ```python
    task.print_summary()
    ```

### Conclusion

Running the simulation on a high-performance machine on the cloud can significantly
reduce the time it takes to complete the simulation. If you need to scale up
or down, you can easly do so by changing the machine to a more or less powerful
one.

Happy simulations!

## Set Commands Manually

In order to give you more flexibility, you can set the commands to run manually,
one by one. This is useful because it gives you more control over the simulation.

If you decide to set your commands manually, here is an example:

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
            on=machine_group)
```

For more details on commands and MPI configuration, refer to the
[Custom Docker Images](https://tutorials.inductiva.ai/simulators/CustomImage.html#command-and-mpiconfig)
documentation.

## What to read next

Try our [Inductiva API](https://console.inductiva.ai/) 
to streamline your workflows and make the most of cloud resources for 
large-scale simulations.

You may also be interested in reading our blog post,
[The 3D Mesh Resolution Threshold - 5k Points is All You Need!](https://inductiva.ai/blog/article/5k-points-is-all-you-need), 
where we explore just how much you can reduce the level of detail in a 
3D object while still maintaining accurate aerodynamic results in a virtual 
wind tunnel built with OpenFOAM.
