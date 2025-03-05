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
12.4 Å³ cell, with MD running for 10 steps.

### Prerequisites  

To follow this tutorial, download the input file for the Fayalite-FIST benchmark from
[here](https://github.com/cp2k/cp2k/blob/master/benchmarks/Fayalite-FIST) and
place it in a folder called `Fayalite-FIST`. Once you have the simulation file, you're
ready to scale your simulations to the Cloud.

### Running Your Simulation

Here is the code required to run a CP2K simulation using the Inductiva API:

```python
"""CP2K Simulation."""
import inductiva

# Instantiate machine group
cloud_machine = inductiva.resources.MachineGroup( 
    provider="GCP",
    machine_type="c3d-standard-16")

# Initialize the Simulator
cp2k = inductiva.simulators.CP2K( 
    version="2025.1")

# Run simulation
task = cp2k.run( 
    input_dir="/Path/to/Fayalite-FIST",
    sim_config_filename="fayalite.inp",
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

We ran this simulation with multiple machines and here are the results.

Starting with a local run on a **Ryzen 7 7700X** with 16 virtual CPUs, the
simulation completed in **1 hour and 3 minutes**. To compare this with a cloud
machine of similar specifications, we used a **c2-standard-16** instance, which
also has 16 vCPUs. However, the cloud machine was slower, taking
**1 hour and 42 minutes**, while costing **0.45 US$**.

| Machine Type            | Virtual CPUs | Time              | Estimated Cost |
|-------------------------|--------------|------------------|---------------|
| **Local Ryzen 7 7700X** | 16           | 1 hour 3 minutes | N/A           |
| **Cloud c2-standard-16** | 16           | 1 hour 42 minutes | 0.45 US$      |

To reduce runtime, we scaled up to a **c2-standard-60** instance with
**60 vCPUs**. This significantly cut the runtime to **42 minutes and 5 seconds**,
but the cost increased to **0.69 US$**.

| Machine Type            | Virtual CPUs | Time              | Estimated Cost |
|------------------------|--------------|------------------|---------------|
| **Cloud c2-standard-16** | 16           | 1 hour 42 minutes | 0.45 US$      |
| **Cloud c2-standard-60** | 60           | 42 minutes 5 seconds | 0.69 US$      |

Instead of just adding more vCPUs, we explored the impact of using newer hardware.
Switching to a **c3d-standard-60** instance—still with 60 vCPUs—improved
performance further, reducing runtime to **32 minutes and 15 seconds**, while
**costing only 0.37 US$**. This highlights that newer hardware can provide better
performance at a lower cost.

| Machine Type            | Virtual CPUs | Time              | Estimated Cost |
|------------------------|--------------|------------------|---------------|
| **Cloud c2-standard-60** | 60           | 42 minutes 5 seconds | 0.69 US$      |
| **Cloud c3d-standard-60** | 60           | 32 minutes 15 seconds | 0.37 US$      |

Finally, we scaled up even further to a **c3d-standard-180** instance with
**180 vCPUs**. The runtime improved slightly to **28 minutes and 12 seconds**,
but the cost rose to **0.96 US$**. While this configuration delivered the best
performance, the additional vCPUs provided diminishing returns in terms of speed
improvement.

| Machine Type            | Virtual CPUs | Time              | Estimated Cost |
|------------------------|--------------|------------------|---------------|
| **Cloud c3d-standard-60** | 60           | 32 minutes 15 seconds | 0.37 US$      |
| **Cloud c3d-standard-180** | 180          | 28 minutes 12 seconds | 0.96 US$      |

### **Final Comparison**
Here's a full overview of all the runs, showing the trade-offs between different setups.

| Machine Type            | Virtual CPUs | Time              | Estimated Cost |
|-------------------------|--------------|------------------|---------------|
| **Local Ryzen 7 7700X** | 16           | 1 hour 3 minutes | N/A           |
| **Cloud c2-standard-16** | 16           | 1 hour 42 minutes | 0.45 US$      |
| **Cloud c2-standard-60** | 60           | 42 minutes 5 seconds | 0.69 US$      |
| **Cloud c3d-standard-60** | 60           | 32 minutes 15 seconds | 0.37 US$      |
| **Cloud c3d-standard-180** | 180          | 28 minutes 12 seconds | 0.96 US$      |

Choosing the right machine for your simulation is crucial to balancing
performance and cost. As seen in the results, blindly increasing vCPUs does not
always lead to proportional speed improvements, and newer hardware can sometimes
provide better efficiency at a lower price. By carefully selecting the right
instance type and size, you can significantly reduce runtime while optimizing
costs, ensuring that your simulations run efficiently without unnecessary expenses.