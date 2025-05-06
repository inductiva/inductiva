# Run the MB9 Micro-benchmark from ExaFOAM
In this tutorial we will show you how to use the Inductiva API to run an advanced OpenFOAM case that requires significant computing power.

## Objective
The goal of this tutorial is to demonstrate how to run the `MB9 Micro-benchmark from ExaFOAM` use case from the CFD Tutorials, available in the [official ExaFOAM documentation](https://exafoam.eu/benchmarks/).

## Prerequisites
Download the required files [here](https://github.com/OpenFOAM/OpenFOAM/tree/25.02/Tutorials/OpenFOAM_CFD/11_2%203D%20Dam%20Break%20with%20Obstacle) and place them in a folder named `highLiftConfiguration`.

The **directory structure** should look like this:
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

## Adjust Simulation Parameters
For a faster simulation, modify the `system/include/caseDefinition` file as follows:
- Set **Time Step** `dt` to `0.00002`
- Set **Start Time** `initTime` to `0.10`
- Set **End Time** `finalTime` to `0.30`
- Set **number of Cores** `nCores` to `180`

Then, you’ll be ready to send your simulation to the Cloud.
 
## Running Your Simulation
Here is the code required to run a OpenFOAM simulation using the Inductiva API:

```python
"""OpenFOAM example"""
import inductiva

# Allocate cloud machine on Google Cloud Platform
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c3d-standard-360",
	spot=True)

# Initialize the Simulator
OpenFOAM = inductiva.simulators.OpenFOAM( \
    version="2412",
	distribution="esi")

# Run simulation
task = OpenFOAM.run(input_dir="/path/to/highLiftConfiguration",
    shell_script="./Allrun",
    on=cloud_machine)

# Wait for the simulation to finish and download the results
task.wait()
cloud_machine.terminate()

task.download_outputs()

task.print_summary()
```

When the simulation is complete, we terminate the machine, download the results and print a summary of the simulation as shown below.

```

```

As you can see in the "In Progress" line, the part of the timeline that represents the actual execution of the simulation, the core computation time 
of this simulation was approximately .

It's that simple! 🚀


