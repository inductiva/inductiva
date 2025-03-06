This guide will walk you through setting up and running CP2K simulations using
the Inductiva API.


# CP2K

[CP2K](https://www.cp2k.org/) is an open-source quantum chemistry and
solid-state physics software designed for atomistic simulations. It supports a
wide range of methods, including Density Functional Theory (DFT), semi-empirical
and classical force fields, making it suitable for molecular
dynamics, electronic structure calculations, and materials science applications.
Its flexibility and extensive features make it a powerful tool for researchers
in computational chemistry, condensed matter physics, and related fields.

## Supported Versions
We currently support the following CP2K version:
- **v2025.1** - Compiled with AVX2.

## Running the CP2K H2O-64 Benchmark

### Objective

This tutorial will show you how to run a CP2K simulation using the H20-64
benchmark available on the [official CP2K website](https://www.cp2k.org/performance#benchmarks).
This use case simulates a system containing 64 water molecules (192 atoms,
512 electrons) in a 12.4 Å³ cell, with MD running for 10 steps.

We will also demonstrate Inductiva's ability to efficiently scale this use case
on a more powerful machine.

### Prerequisites  

To follow this tutorial, download the input file for the H2O-64 benchmark from
[here](https://github.com/cp2k/cp2k/blob/master/benchmarks/QS/H2O-64.inp) and
place it in a folder called `H2O-64`. Once you have the simulation file, you're
ready to scale your simulations to the Cloud.


### Running Your Simulation

Here is the code required to run a CP2K simulation using the Inductiva API.

We will be running this simulation on a 16 vCPU virtual machine supported by a
4th generation AMD EPYC™ (Genoa) processor.

```python
"""CP2K Simulation."""
import inductiva

# Instantiate machine group
cloud_machine = inductiva.resources.MachineGroup( 
    provider="GCP",
    machine_type="c3d-highcpu-16")

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
inductiva tasks info 7y33ri1umg0tsppabxf5jska4

Task status: Success

Timeline:
	Waiting for Input         at 05/03, 16:04:17      1.525 s
	In Queue                  at 05/03, 16:04:19      21.308 s
	Preparing to Compute      at 05/03, 16:04:40      8.307 s
	In Progress               at 05/03, 16:04:48      103.301 s
		└> 103.177 s       /opt/openmpi/5.0.6/bin/mpirun --use-hwthread-cpus --np 16 cp2k.psmp H2O-64.inp
	Finalizing                at 05/03, 16:06:32      0.453 s
	Success                   at 05/03, 16:06:32      

Data:
	Size of zipped output:    87.62 KB
	Size of unzipped output:  296.14 KB
	Number of output files:   6

Estimated computation cost (US$): 0.0049 US$

Go to https://console.inductiva.ai/tasks/7y33ri1umg0tsppabxf5jska4 for more details.
```

The core computation time for this simulation was approximately **1 minute and 43 seconds**
(103 seconds), as shown in the `In Progress` line. This represents the
actual execution time of the CP2K benchmark on a 16 virtual CPU machine.

For comparison, the same simulation takes **1 minute and 15 seconds** on a similar
local machine with a 16 virtual CPUs (Ryzen 7 7700X). This performance
difference is expected, as cloud CPUs typically have lower clock speeds compared to
high-performance desktop processors, prioritizing energy efficiency and density
over raw speed.

However, increasing the number of vCPUs on this cloud machine can improve this
result.

### Scaling Up Your Simulation  

With the Inductiva API, scaling up your CP2K simulation is as simple as changing
two parameters:

1. Modify the `machine_type` to a more powerful machine with more vCPUs.
2. Adjust the `n_vcpus` accordingly to maximize parallel processing efficiency.

We tested this simulation across multiple machines to analyze how performance
and cost scale with increasing computational resources.  

We began with a local run on a **Ryzen 7 7700X** with **16 vCPUs**, completing
the simulation in **1 minute and 15 seconds**. To compare this with a
cloud-based machine of similar specifications, we used a **c3d-highcpu-16**
machine, which also has **16 vCPUs**. As expected, the cloud machine was a bit slower,
taking **1 minute and 43 seconds**, with a cost of **0.0049 US$**.

| Machine Type            | Virtual CPUs | Time              | Estimated Cost |
|-------------------------|--------------|------------------|---------------|
| **Local Ryzen 7 7700X** | 16           | 1 minute and 15 seconds | N/A           |
| **Cloud c3d-highcpu-16** | 16           | 1 minute and 43 seconds | 0.0049 US$      |

To improve performance, we scaled up to a **c3d-highcpu-60** machine with
**60 vCPUs**. This significantly reduced the runtime to **38 seconds**, with the
cost increasing slightly to **0.0077 US$**.  

| Machine Type            | Virtual CPUs | Time              | Estimated Cost |
|------------------------|--------------|------------------|---------------|
| **Cloud c3d-highcpu-16** | 16           | 1 minute and 43 seconds | 0.0049 US$      |
| **Cloud c3d-highcpu-60** | 60           | 38 seconds | 0.0077 US$      |

### **Final Comparison**  

| Machine Type            | Virtual CPUs | Time              | Estimated Cost |
|-------------------------|--------------|------------------|---------------|
| **Local Ryzen 7 7700X** | 16           | 1 minute and 15 seconds | N/A           |
| **Cloud c3d-highcpu-16** | 16           | 1 minute and 43 seconds | 0.0049 US$      |
| **Cloud c3d-highcpu-60** | 60           | 38 seconds | 0.0077 US$      | 

By leveraging the Inductiva API, you can efficiently scale your CP2K simulations
to meet your computational needs. Try different machine configurations and
optimize your workflow for faster and more cost-effective results!