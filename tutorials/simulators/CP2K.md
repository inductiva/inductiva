This guide will walk you through running your first CP2K simulations using the Inductiva API.


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

This tutorial will demonstrate how to run the H2O-64 benchmark, which simulates
a system that consists of 64 water molecules in a 12.4 Å³ cell, with MD running
for 10 steps.

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

For comparison, this same simulation takes **1 minute and 15 seconds** on a similar
local machine with a 16 virtual CPUs (Ryzen 7 7700X). This performance
difference is expected, as cloud CPUs typically have lower clock speeds compared to
high-performance desktop processors, prioritizing energy efficiency and density
over raw speed.

The good thing about Inductiva is that you only need to change two lines of code
to speed up your simulation.

### Scaling Up Your Simulation  

Scaling up your CP2K simulation is as simple as changing two parameters:

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

### How much is too much?

To further explore cloud scaling, we tested two additional machines:
**c3d-highcpu-180** and **c3d-highcpu-360**. However, for a small simulation
like this, scaling up to such high vCPU counts does not necessarily yield better
performance. The results show a slight improvement for the **c3d-highcpu-180** machine,
with a  time of **33 seconds** and cost of **0.020 US$**, and a slow down for
the **c3d-highcpu-360** machine, with a time of **34 seconds** and cost of **0.045 US$**.

| Machine Type            | Virtual CPUs | Time              | Estimated Cost |
|------------------------|--------------|------------------|---------------|
| **Cloud c3d-highcpu-180** | 180           | 33 seconds | 0.020 US$      |
| **Cloud c3d-highcpu-360** | 360           | 34 seconds | 0.045 US$      |

### **Final Comparison**  

| Machine Type            | Virtual CPUs | Time              | Estimated Cost |
|-------------------------|--------------|------------------|---------------|
| **Local Ryzen 7 7700X** | 16           | 1 minute and 15 seconds | N/A           |
| **Cloud c3d-highcpu-16** | 16           | 1 minute and 43 seconds | 0.0049 US$      |
| **Cloud c3d-highcpu-60** | 60           | 38 seconds | 0.0077 US$      |
| **Cloud c3d-highcpu-180** | 180           | 33 seconds | 0.020 US$      |
| **Cloud c3d-highcpu-360** | 360           | 34 seconds | 0.045 US$      |

### **Key Takeaway**  

Selecting the right machine for your simulation is essential for balancing
performance and cost. As the results show, simply increasing vCPUs does not
always lead to faster execution, and in some cases, it may even slow down the
simulation. Ultimately, the best choice depends on your priorities—whether you
aim for the fastest possible runtime, the lowest cost, or a balanced approach
that delivers good performance at a reasonable price.
