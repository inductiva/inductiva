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

To follow this tutorial, download the input file for the H2O-64 benchmark from
[here](https://github.com/cp2k/cp2k/blob/master/benchmarks/QS/H2O-64.inp) and
place it in a folder called `H2O-64`. Once you have the simulation file, you're
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
inductiva tasks info ngjrax2xsak3yaonirxj1pdcc

Task status: Success

Timeline:
	Waiting for Input         at 05/03, 10:28:59      1.193 s
	In Queue                  at 05/03, 10:29:00      12.325 s
	Preparing to Compute      at 05/03, 10:29:13      6.645 s
	In Progress               at 05/03, 10:29:19      103.302 s
		└> 103.17 s        /opt/openmpi/4.1.6/bin/mpirun --use-hwthread-cpus --np 16 cp2k.psmp H2O-64.inp
	Finalizing                at 05/03, 10:31:03      0.447 s
	Success                   at 05/03, 10:31:03      

Data:
	Size of zipped output:    87.65 KB
	Size of unzipped output:  296.30 KB
	Number of output files:   6

Estimated computation cost (US$): 0.0058 US$

Go to https://console.inductiva.ai/tasks/ngjrax2xsak3yaonirxj1pdcc for more details.
```

The core computation time for this simulation was approximately 1 minute and 43
seconds (103 seconds), as shown in the `In Progress` line. This represents the
actual execution time of the CP2K benchmark on a 16 virtual CPU machine.

### Scaling Up Your Simulation  

Scaling up your CP2K simulation is as simple as changing two parameters:

1. Modify the `machine_type` to a more powerful machine with more vCPUs.
2. Adjust the `n_vcpus` accordingly to maximize parallel processing efficiency.

We tested this simulation across multiple machines to analyze how performance
and cost scale with increasing computational resources.  

We began with a local run on a **Ryzen 7 7700X** with **16 vCPUs**, completing
the simulation in **1 minute and 15 seconds**. To compare this with a
cloud-based machine of similar specifications, we used a **c3d-standard-16**
machine, which also has **16 vCPUs**. However, the cloud machine was slower,
taking **1 minute and 50 seconds**, with a cost of **0.0058 US$**. This
performance difference is expected, as cloud CPUs typically have lower clock
speeds compared to high-performance desktop processors, prioritizing energy
efficiency and density over raw speed.

| Machine Type            | Virtual CPUs | Time              | Estimated Cost |
|-------------------------|--------------|------------------|---------------|
| **Local Ryzen 7 7700X** | 16           | 1 minute and 15 seconds | N/A           |
| **Cloud c3d-standard-16** | 16           | 1 minute and 50 seconds | 0.0058 US$      |

To improve performance, we scaled up to a **c3d-standard-60** machine with
**60 vCPUs**. This significantly reduced the runtime to **51 seconds**, with the
cost increasing slightly to **0.0098 US$**.  

| Machine Type            | Virtual CPUs | Time              | Estimated Cost |
|------------------------|--------------|------------------|---------------|
| **Cloud c3d-standard-16** | 16           | 1 minute and 50 seconds | 0.0058 US$      |
| **Cloud c3d-standard-60** | 60           | 51 seconds | 0.0098 US$      |

To further explore cloud scaling, we tested two additional machines:
**c3d-standard-180** and **c3d-standard-360**. However, for a small simulation
like this, scaling up to such high vCPU counts does not necessarily yield better
performance. The results show worst results, with runtimes increasing to
**1 minute and 5 seconds** and **1 minute and 48 seconds**, while costs surged
to **0.037 US$** and **0.12 US$**, respectively.  

| Machine Type            | Virtual CPUs | Time              | Estimated Cost |
|------------------------|--------------|------------------|---------------|
| **Cloud c3d-standard-180** | 180           | 1 minute and 5 seconds | 0.037 US$      |
| **Cloud c3d-standard-360** | 360           | 1 minute and 48 seconds | 0.12 US$      |

### **Final Comparison**  

| Machine Type            | Virtual CPUs | Time              | Estimated Cost |
|-------------------------|--------------|------------------|---------------|
| **Local Ryzen 7 7700X** | 16           | 1 minute and 15 seconds | N/A           |
| **Cloud c3d-standard-16** | 16           | 1 minute and 50 seconds | 0.0058 US$      |
| **Cloud c3d-standard-60** | 60           | 51 seconds | 0.0098 US$      |
| **Cloud c3d-standard-180** | 180           | 1 minute and 5 seconds | 0.037 US$      |
| **Cloud c3d-standard-360** | 360           | 1 minute and 48 seconds | 0.12 US$      |

### **Key Takeaway**  

Selecting the right machine for your simulation is essential for balancing
performance and cost. As the results show, simply increasing vCPUs does not
always lead to faster execution, and in some cases, it may even slow down the
simulation. Ultimately, the best choice depends on your priorities—whether you
aim for the fastest possible runtime, the lowest cost, or a balanced approach
that delivers good performance at a reasonable price.
