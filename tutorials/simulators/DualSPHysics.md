In this guide, we will walk you through setting up and running DualSPHysics, 
a Smoothed-Particle Hydrodynamics (SPH) simulator, available as one of 
the built-in tools via the Inductiva API. 

We will cover:

- Setting up DualSPHysics for use with our API.
- Example code to help you get started with simulations.
- An advanced Turbine example to show how to execute commands through the 
Inductiva API.

# DualSPHysics

DualSPHysics is a Smoothed-Particle Hydrodynamics (SPH) simulator. The simulator 
is usually configured by a single file with the extension `.xml`. This file
contains all the information about the simulation, including the geometry, 
the physical properties of the fluids, the boundary conditions, the numerical
parameters, and the output files. Sometimes the configuration can also use extra
geometry files. 

Running your DualSPHysics simulation workflows using Inductiva is very similar
to running them on your local machine, but instead of calling your DualSPHysics
shell script directly, you will have to pass it to the Inductiva API via the
run() method to be executed on a remote resource. So, if you already have a
functioning shell script that orchestrates the entire DualSPHysics simulation on
your local machine, you are mostly almost ready to run it via Inductiva.

There are, however, some minor adaptations that may have to be done to your
orchestration script to take into account the difference in the environment.
For example, you may have to change path-related variables in your script
because the location of the DualSPHysics binaries in our infrastructure may be
different from that of your local setup.

Also, interactive commands that wait for your keyboard input may have to be
removed or set to run on default parameters.

## Technical details

This section will focus on the implementation details of our compilation of 
DualSPHysics.

Lets start with the available versions of DualSPHysics. For that you can check
our dockerhub page [here](https://hub.docker.com/r/inductiva/kutu/tags?name=dualsphysics).

**Location of binaries**: All binaries related to DualSPHysics are located in
`/DualSPHysics_v5.2/bin/linux/`. This path was added to the environment variable
`PATH`, so you can call the commands directly without the need to specify the full
path if you wish so. We also created a list of simbolic links to allow you to
call some binaries without the need to know the full name of the binary. For
example, you can call `dualsphysics` instead of `DualSPHysics5.2CPU_linux64`.
This renaming follows a pattern that is easy to understand. We removed the
`_linux64` suffix and make the name lowercase. So, you can either use the names
you are used to (`DualSPHysics5.2CPU_linux64`) or use the simplified names that
will abstract achitecture and simulator version (`dualsphysics`).

Below, we have a concrete example that will let you better understand the
changes you may potentially have to do.

For an extensive list of commands, please refer to the DualSPHysics 
[documentation](https://dual.sphysics.org/). You can pass the API commands in
lowercase, and we will handle the rest for you!

## Example code

In this example, we run a classical CFD case of a flow over a cylinder. 

```{literalinclude} ../../examples/dualsphysics/dualsphysics.py
:language: python
```

## Advanced Tutorial: Running `examples/chrono/09_Turbine`

### Objective
This tutorial demonstrates how to run the advanced `09_Turbine` example included
in the DualSPHysics distribution using the Inductiva API.

### Prerequisites

1. **Download Input Files**: Download DualSPHysics [package](https://dual.sphysics.org/downloads/)
and see if you can find the example `examples/chrono/09_Turbine`. Navigate to 
`examples/chrono`. We are going to work from that directory and write our
Inductiva python script there.

2. **Update Simulation Script**: Update the simulation script of the example `09_Turbine`.

### Overview
The first step involves making changes to the `xCaseTurbine_linux64_CPU.sh` script for
execution within our environment. This step is straightforward and requires only
minor modifications to the original script. The simulation generates output files for 
visualization in ParaView.

Here is the overview of the code for this simulation:

```python
"""DualSPHysics example."""
import inductiva

# Instantiate machine group
machine_group = inductiva.resources.MachineGroup(
    machine_type="n2d-highcpu-64",
    spot=True,
    data_disk_gb=200)
machine_group.start()

# Download the configuration files into a folder
input_dir = "/Path/to/09_Turbine"


# Initialize the Simulator
dualsphysics = inductiva.simulators.DualSPHysics()

# Run simulation with config files in the input directory
task = dualsphysics.run(input_dir=input_dir,
                        shell_script="xCaseTurbine_linux64_CPU.sh",
                        on=machine_group)

task.wait()
task.download_outputs()

machine_group.terminate()
```

### Step 1: Adjust Simulation Script

Before running the simulation, we need to adjust the simulation script located
at `examples/chrono/09_Turbine/xCaseTurbine_linux64_CPU.sh`.

1. **Update the `dirbin` Variable:**
   Modify the `xCaseTurbine_linux64_CPU.sh` script to point to the correct binaries directory:
   ```bash
   export dirbin=/DualSPHysics_v5.2/bin/linux/
   ```

2. **Remove User Input Prompt:**
   Remove the last line on the script. This line waits for user input and will
   prevent the script from running in an automated environment:
   
   ```bash
   read -n1 -r -p "Press any key to continue..." key
   ```

These modifications prepare the script for automated execution.

### Step 2: Running the Simulation

#### a. Configure and Start Machine

In order to chose the right machine for your simulation, you need to be aware
of our computational infrastructure. We sujest you to read the documentation
[here](https://docs.inductiva.ai/en/latest/intro_to_api/computational-infrastructure.html).

For this simulation we decided to go with a `n2d-highcpu-64` machine. This
machine has 64 virtual CPUs and a 20 GB data disk. We also decided to use a spot
machine to reduce the cost of the simulation.

```python
import inductiva

machine_group = inductiva.resources.MachineGroup(
      machine_type="n2d-highcpu-64",
      spot=True,
      data_disk_gb=20)

machine_group.start()
```

#### b. Simulation inputs

1. **Specify Simulation Directory**:
   - Set the `input_dir` parameter to point to the example folder `09_Turbine`.
   ```python
    input_dir = "/Path/to/09_Turbine"
   ```

#### c. Run your simulation

1. **Pick the simulator**:
   - Initialize the DualSPHysics simulator:
   ```python
   dualsphysics = inductiva.simulators.DualSPHysics()
   ```

2. **Run the simulation**:
   Execute the simulation with the prepared script and input files:
   ```python
   task = dualsphysics.run(
       input_dir=input_dir,
       shell_script="xCaseTurbine_linux64_CPU.sh",
       on=machine_group)
   ```

2. **Wait**:
   Wait for the simulation to complete:
   ```python
   task.wait()
   ```
   > Note: The `wait()` method will block the execution of the script until the
   simulation is completed. This is useful for scripts that need to perform
   additional actions after the simulation completes.

3. **Terminate Machine**:
   After the simulation completes, terminate the machine group:
   ```python
   machine_group.terminate()
   ```

4. **Check your simulation summary**:
   View the task summary with:
   ```python
   task.print_summary()
   ```

   Wich will output something like:
   ```bash
    ■ Tier: Power-User

    ■ Credits: 1000.00 US$

    ■ Global User quotas
                                                                    CURRENT USAGE     MAX ALLOWED
    Maximum simultaneous instances                                  0 instance        100 instance
    Maximum price per hour across all instances                     0 USD             270 USD
    Maximum tasks per week                                          3 task            N/A
    Maximum number of VCPUs                                         0 vcpu            1000 vcpu
    Maximum time a machine group can stay idle before termination   N/A               120 minute

    ■ Instance User quotas
                                                                                              MAX ALLOWED
    Maximum disk size                                                                        2000 GB
    Maximum time a machine group can stay up before automatic termination                    48 hour
    Maximum amount of RAM per VCPU                                                           6 GB

    ■ Registering MachineGroup configurations:
      · Name:                       api-57yp64vdyz50a8fx0ttuz2n90
      · Machine Type:               n2d-highcpu-64
      · Data disk size:             20 GB
      · Maximum idle time:          30 minutes
      · Auto terminate timestamp:   2024/08/20 20:01:22
      · Number of machines:         1
      · Spot:                       True
      · Estimated cloud cost of machine group: 0.670 $/h
      · You are spending 3.3x less by using spot machines.

    Starting MachineGroup(name="api-57yp64vdyz50a8fx0ttuz2n90"). This may take a few minutes.
    Note that stopping this local process will not interrupt the creation of the machine group. Please wait...
    Machine Group api-57yp64vdyz50a8fx0ttuz2n90 with n2d-highcpu-64 machines successfully started in 0:00:29.

    The machine group is using the following quotas:

                                                  USED BY RESOURCE     CURRENT USAGE     MAX ALLOWED
    Maximum number of VCPUs                       64                   64                1000
    Maximum simultaneous instances                1                    1                 100
    Maximum price per hour across all instances   0.6727               0.6727            270

    ■ Using production image of DualSPHysics version 5.2.1

    ■ Task Information:
      · ID:                    u8v7p1v7wfyvvkyc0iq0s632k
      · Simulator:             DualSPHysics
      · Version:               5.2.1
      · Image:                 docker://inductiva/kutu:dualsphysics_v5.2.1
      · Local input directory: examples/chrono/09_Turbine
      · Submitting to the following computational resources:
        · Machine Group api-57yp64vdyz50a8fx0ttuz2n90 with n2d-highcpu-64 machines

    Preparing upload of the local input directory examples/chrono/09_Turbine (2.25 MB).
    Input archive size: 1.09 MB
    Uploading input archive...
    100%|██████████████████████████████████████████████████████████████████████████████| 1.09M/1.09M [00:01<00:00, 715kB/s]
    Local input directory successfully uploaded.

    ■ Task u8v7p1v7wfyvvkyc0iq0s632k submitted to the queue of the Machine Group api-57yp64vdyz50a8fx0ttuz2n90 with n2d-highcpu-64 machines.
    Number of tasks ahead in the queue: 0
    · Consider tracking the status of the task via CLI:
      inductiva tasks list --id u8v7p1v7wfyvvkyc0iq0s632k
    · Or, tracking the logs of the task via CLI:
      inductiva logs u8v7p1v7wfyvvkyc0iq0s632k
    · You can also get more information about the task via the CLI command:
      inductiva tasks info u8v7p1v7wfyvvkyc0iq0s632k

    Waiting for task u8v7p1v7wfyvvkyc0iq0s632k to complete...
    Go to https://console.inductiva.ai/tasks/u8v7p1v7wfyvvkyc0iq0s632k for more details.
    ■ Task u8v7p1v7wfyvvkyc0iq0s632k successfully queued and waiting to be picked-up for execution...
    The task u8v7p1v7wfyvvkyc0iq0s632k is about to start.
    ■ Task u8v7p1v7wfyvvkyc0iq0s632k has started and is now running remotely.
    ■ Task u8v7p1v7wfyvvkyc0iq0s632k completed successfully.
    Downloading stdout and stderr files to u8v7p1v7wfyvvkyc0iq0s632k...
    Partial download completed to u8v7p1v7wfyvvkyc0iq0s632k.
    Successfully requested termination of MachineGroup(name="api-57yp64vdyz50a8fx0ttuz2n90").
    Termination of the machine group freed the following quotas:

                                                  FREED BY RESOURCE     CURRENT USAGE     MAX ALLOWED
    Maximum number of VCPUs                       64                    0                 1000
    Maximum simultaneous instances                1                     0                 100
    Maximum price per hour across all instances   0.6727                0                 270


    Task status: success
    Wall clock time:  0:09:35
    Time breakdown:
      Input upload:              1.81 s
      Time in queue:             10.71 s
      Container image download:  1.23 s
      Input download:            0.09 s
      Input decompression:       0.01 s
      Computation:               0:06:20
      Output upload:             0:03:01
    Data:
      Size of zipped output:    3.52 GB
      Size of unzipped output:  5.35 GB
      Number of output files:   2541

   ```

That's it!

We can now donwload the results to our local machine using Inductiva's CLI:
```bash
inductiva tasks download u8v7p1v7wfyvvkyc0iq0s632k
```
Downloading and decompressing data will take a few minutes (depending on your
network connection):
```bash
Downloading simulation outputs to inductiva_output/u8v7p1v7wfyvvkyc0iq0s632k/output.zip...
100%|█████████████████████████████████████████████████████████████████████████████| 3.52G/3.52G [04:43<00:00, 12.4MB/s]
Uncompressing the outputs to u8v7p1v7wfyvvkyc0iq0s632k...
```

As usual, the results are placed in the `inductiva_output` folder, within a 
subfolder named after the task. Earlier, we set a variable for the internal
directory where all outputs would be placed (`dirout`), which was instantiated 
as `CaseTurbine_out`. Let’s check its contents:

```bash
ls -las inductiva_output/u8v7p1v7wfyvvkyc0iq0s632k/CaseTurbine_out
total 36080
   0 drwxr-xr-x    22 lsarmento  staff      704 19 Aug 12:09 .
   0 drwxr-xr-x    14 lsarmento  staff      448 19 Aug 12:09 ..
9888 -rw-r--r--     1 lsarmento  staff  5058947 19 Aug 12:07 CaseTurbine.bi4
  16 -rw-r--r--     1 lsarmento  staff     4523 19 Aug 12:07 CaseTurbine.out
  24 -rw-r--r--     1 lsarmento  staff    10830 19 Aug 12:07 CaseTurbine.xml
9432 -rw-r--r--     1 lsarmento  staff  4827935 19 Aug 12:07 CaseTurbine_All.vtk
1280 -rw-r--r--     1 lsarmento  staff   653967 19 Aug 12:07 CaseTurbine_Bound.vtk
8160 -rw-r--r--     1 lsarmento  staff  4174240 19 Aug 12:07 CaseTurbine_Fluid.vtk
 752 -rw-r--r--     1 lsarmento  staff   382901 19 Aug 12:07 CaseTurbine_MkCells.vtk
2488 -rw-r--r--     1 lsarmento  staff  1272363 19 Aug 12:07 CaseTurbine__Actual.vtk
   8 -rw-r--r--     1 lsarmento  staff      583 19 Aug 12:07 CaseTurbine_dbg-fillbox.vtk
   8 -rw-r--r--     1 lsarmento  staff     1947 19 Aug 12:07 CfgChrono_Scheme.vtk
   8 -rw-r--r--     1 lsarmento  staff      854 19 Aug 12:07 CfgInit_Domain.vtk
  16 -rw-r--r--     1 lsarmento  staff     4155 19 Aug 12:07 CfgInit_MapCells.vtk
   8 -rw-r--r--     1 lsarmento  staff     2415 19 Aug 12:07 Floating_Materials.xml
3880 -rw-r--r--     1 lsarmento  staff  1985884 19 Aug 12:07 Rotor.stl
   8 -rw-r--r--     1 lsarmento  staff      909 19 Aug 12:07 Run.csv
 104 -rw-r--r--     1 lsarmento  staff    49764 19 Aug 12:07 Run.out
   0 drwxr-xr-x   504 lsarmento  staff    16128 19 Aug 12:09 boundary
   0 drwxr-xr-x   508 lsarmento  staff    16256 19 Aug 12:08 data
   0 drwxr-xr-x  1006 lsarmento  staff    32192 19 Aug 12:08 particles
   0 drwxr-xr-x   503 lsarmento  staff    16096 19 Aug 12:09 surface
```
Data for visualization is placed inside the directories `boundary`, `particles`
and `surface`. This data can be loaded in ParaView and rendered in a movie as 
the one seen above.

