This guide will walk you through setting up and running CP2K simulations using the Inductiva API.


# CP2K

[CP2K](https://www.cp2k.org/) is an open-source quantum chemistry and
solid-state physics software designed for atomistic simulations. It supports a
wide range of methods, including density functional theory (DFT), semi-empirical
approaches, and classical force fields, making it suitable for molecular
dynamics, electronic structure calculations, and materials science applications.
Its flexibility and extensive features make it a powerful tool for researchers
in computational chemistry, condensed matter physics, and related fields.

## Supported Versions
We currently support the following CP2K version:
- **v2025.1** - Compiled with AVX2.

## Running the CP2K H2O-64 Benchmark

### Objective

This tutorial will demonstrate how to run a CP2K simulation using the H2O-64
benchmark, which simulates a system that consists of 64 water molecules in a
12.4 Å³ cell, with MD running for 1000 steps.

### Prerequisites  

To follow this tutorial, download the input file for the H2O-64 benchmark from
[here](https://github.com/cp2k/cp2k/blob/master/benchmarks/QS/H2O-64.inp) and
place it in a folder called `H2O-64`. Once you have the simulation file, you're
ready to scale your simulations to the Cloud.

### Increasing the simulation time

The H2O-64 benchmark only runs for 10 steps as configured, which corresponds to
a simulated time of 5 femtoseconds. For the purposes of this tutorial, we
decided to increase the simulation time by a factor of 100 to 500 femtoseconds.
This requires increasing the number of MD steps from 10 to 1000 in the CP2K
input file.

To do this, open the file `H2O-64.inp` and change `STEPS 10` to `STEPS 1000`.

### Running Your Simulation

Here is the code required to run a CP2K simulation using the Inductiva API:

```python
"""CP2K Simulation."""
import inductiva

# Instantiate machine group
cloud_machine = inductiva.resources.MachineGroup( 
    provider="GCP",
    machine_type="c2-standard-16")

# Initialize the Simulator
cp2k = inductiva.simulators.CP2K( 
    version="2025.1")

# Run simulation
task = cp2k.run( 
    input_dir="/Path/to/H2O-64",
    sim_config_filename="H2O-64.inp",
    n_vcpus=16,
    use_hwthread=True,
    on=cloud_machine)

task.wait()
cloud_machine.terminate()

task.download_outputs()

task.print_summary()
```

To adapt this script for other CP2K simulations, replace `input_dir` with the
path to your CP2K input files and set the `sim_config_filename` accordingly.

Once the simulation is complete, we terminate the machine, download the results,
and print a summary of the simulation as shown below.

```
inductiva tasks info 73bneshwui82w5entxwxjk9j7

Task status: Success

Timeline:
	Waiting for Input         at 04/03, 09:21:50      1.032 s
	In Queue                  at 04/03, 09:21:51      34.652 s
	Preparing to Compute      at 04/03, 09:22:26      10.603 s
	In Progress               at 04/03, 09:22:36      6134.309 s
		└> 6134.176 s      /opt/openmpi/4.1.6/bin/mpirun --use-hwthread-cpus --np 16 cp2k.psmp H2O-64.inp
	Finalizing                at 04/03, 11:04:51      0.956 s
	Success                   at 04/03, 11:04:51      

Data:
	Size of zipped output:    5.33 MB
	Size of unzipped output:  16.61 MB
	Number of output files:   9

Estimated computation cost (US$): 0.46 US$

Go to https://console.inductiva.ai/tasks/73bneshwui82w5entxwxjk9j7 for more details.
```

The core computation time for this simulation was approximately 1 hour and 42
minutes (6134 seconds), as shown in the `In Progress` line. This represents the
actual execution time of the CP2K benchmark on a 16 virtual CPU machine.

### Scaling Up Your Simulation  

Scaling up your CP2K simulation is as simple as changing two parameters:

1. Modify the `machine_type` to a more powerful machine with more vCPUs.
2. Adjust the `n_vcpus` accordingly to maximize parallel processing efficiency.

Here are the results of running the H2O-64 benchmark on different machines:

|  Machine Type  | Virtual CPUs |     Time     | Estimated Cost |
|:--------------:|:------------:|:------------:|:--------------:|
|  Local Ryzen 7 7700X |      16      | 1 hour and 3 minutes       | N/A       |
|  Cloud c2-standard-16      |      16      | 1 hour and 42 minutes      | 0.45 US$   |
|  Cloud c2-standard-60      |      60      | 42 minutes and 5 seconds   | 0.69 US$   |
|  Cloud c3d-standard-60      |      60      | 32 minutes and 15 seconds   | 0.37 US$   |
|  Cloud c3d-standard-180      |      180      | 28 minutes and 12 seconds   | 0.96 US$   |

Running the CP2K simulation on a local Ryzen 7 7700X with 16 cores took 1 hour
and 3 minutes as the baseline. When moving to a similar cloud machine
(c2-standard-16) with 16 vCPUs, the simulation took longer—1 hour and 42
minutes—likely due to lower clock speeds, but at a low cost of 0.45 US$.  

Scaling up to a more powerful cloud machine (c3d-standard-60) with 60 vCPUs
significantly reduced the simulation time to 32 minutes and 15 seconds while
decreasing the cost to 0.37 US$.

This highlights the importance of choosing the right machine for your simulation:
we more than halved the computation time while also decreasing the cost of the simulation.
