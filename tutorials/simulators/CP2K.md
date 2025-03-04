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
- **v2025.1**

## Running the CP2K H2O-64 Benchmark

### Objective

This tutorial will demonstrate how to run a CP2K simulation using the H2O-64
benchmark, which simulates a system that consists of 64 water molecules in a
12.4 Å³ cell, with MD running for 1000 steps.

### Prerequisites  

To follow this tutorial, download the input file for the H2O-64 benchmark from
[here](https://github.com/cp2k/cp2k/blob/master/benchmarks/QS/H2O-64.inp) and
place it in a folder called `H2O-64`. Once you have the simulation file, you're
ready to scale your simulations to the cloud.

### Increasing the simulation time

The H2O-64 benchmark, as configured, runs for only 10 steps, which corresponds
to a simulated time of 5 femtoseconds. To obtain more relevant results, we
increase the simulation time by 100x, extending it to 500 femtoseconds. This
requires increasing the number of MD steps from 10 to 1000 in the CP2K input file.

To do that open the `H2O-64.inp` file and change `STEPS 10` to `STEPS 1000`.

### Running Your Simulation

Here is the code required to run a CP2K simulation using the Inductiva API:

```python
"""CP2K Simulation."""
import inductiva

# Instantiate machine group
cloud_machine = inductiva.resources.MachineGroup( 
    provider="GCP",
    machine_type="c2-standard-4")

# Initialize the Simulator
cp2k = inductiva.simulators.CP2K( 
    version="2025.1")

# Run simulation
task = cp2k.run( 
    input_dir="/Path/to/H2O-64",
    sim_config_filename="H2O-64.inp",
    n_vcpus=4,
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
inductiva tasks info 6qcy46uvjdoqyysx8d0zbjv5w

Task status: Success

Timeline:
	Waiting for Input         at 03/03, 20:00:16      1.057 s
	In Queue                  at 03/03, 20:00:17      12.13 s
	Preparing to Compute      at 03/03, 20:00:29      8.826 s
	In Progress               at 03/03, 20:00:38      22471.079 s
		└> 22470.94 s      /opt/openmpi/4.1.6/bin/mpirun --use-hwthread-cpus --np 4 cp2k.psmp H2O-64.inp
	Finalizing                at 04/03, 02:15:09      0.766 s
	Success                   at 04/03, 02:15:09      

Data:
	Size of zipped output:    5.33 MB
	Size of unzipped output:  16.58 MB
	Number of output files:   9

Estimated computation cost (US$): 0.44 US$

Go to https://console.inductiva.ai/tasks/6qcy46uvjdoqyysx8d0zbjv5w for more details.
```

The core computation time for this simulation was approximately 6 hours and 14
minutes (22471 seconds), as shown in the `In Progress` line. This represents the
actual execution time of the CP2K benchmark on a 4 virtual CPU machine.

### Scaling Up Your Simulation  

Scaling up your CP2K simulation is as simple as changing two parameters:

1. Modify the `machine_type` to a more powerful machine with more vCPUs.
2. Adjust the `n_vcpus` accordingly to maximize parallel processing efficiency.

Here are the results of running the H2O-64 benchmark on different machines:

|  Machine Type  | Virtual CPUs |     Time     | Estimated Cost |
|:--------------:|:------------:|:------------:|:--------------:|
|  Local Ryzen 7 7700X |      16      | ...   | N/A       |
|  c2-standard-4 |      4      | 6 hours and 14 minutes   | 0.43 US$       |
|  c2-standard-16 |      16      | N/A  | N/A       |
|  c2-standard-60 |      60      | 42 minutes and 5 seconds   | 0.69 US$   |

We also ran the present simulation using GPU capable machines so see what kind
of speed ups we get compared to running on the CPU only.

|  Machine Type  | Virtual CPUs | GPU          |     Time     | Estimated Cost |
|:--------------:|:------------:|:------------:|:------------:|:--------------:|
|  Local RTX 4070|      16      | 1x RTX 4070 | 6 hours and 14 minutes   | N/A |
|  g2-standard-4 |      4       | 1x NVIDIA L4 | 6 hours and 14 minutes   | 0.43 US$       |



CP2K simulations can be computationally intensive, but with the right hardware,
you can significantly speed up your simulations and reduce overall cost.

