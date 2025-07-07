# Run Your First Simulation
This tutorial will show you how to run CP2K simulations using the Inductiva API. 

This tutorial will cover the `H20-64 benchmark`, available in the [official CP2K website](https://www.cp2k.org/performance#benchmarks), to help you get started with simulations. 

The use case simulates a system containing 64 water molecules (192 atoms, 512 electrons) in a 12.4 Å³ cell, with MD running for 10 steps.

We will also demonstrate Inductiva’s ability to efficiently scale this use case on a more powerful machine.

<p align="center"><img src="./_static/h2o-64.gif" alt="CP2K simulation visualization" width="700"></p>

## Prerequisites
Download the required files [here](https://github.com/cp2k/cp2k/blob/v2025.1/benchmarks/QS/H2O-64.inp) and place them in a folder called `H2O-64`. Then, you’ll be ready to send your simulation to the Cloud.

## Running an CP2K Simulation
Here is the code required to run OpenSees simulation using the Inductiva API:

```python
"""CP2K Simulation."""
import inductiva

# Allocate a machine on Google Cloud Platform
cloud_machine = inductiva.resources.MachineGroup( 
    provider="GCP",
    machine_type="c2d-highcpu-16",
    spot=True)

# Initialize the Simulator
cp2k = inductiva.simulators.CP2K( 
    version="2025.1")

# Run simulation
task = cp2k.run( 
    input_dir="/Path/to/H2O-64", 
    sim_config_filename="H2O-64.inp",
    on=cloud_machine)

# Wait for the simulation to finish and download the results
task.wait()
cloud_machine.terminate()

task.download_outputs()

task.print_summary()
```

In this basic example, we're using a cloud machine (`c2d-highcpu-16`) equipped with 16 virtual CPUs. 
For larger or more compute-intensive simulations, consider adjusting the `machine_type` parameter to select 
a machine with more virtual CPUs or one equipped with GPUs. You can explore the full range of available 
machines [here](https://console.inductiva.ai/machine-groups/instance-types).

> **Note**: Setting `spot=True` enables the use of spot machines, which are available at substantial discounts. 
> However, your simulation may be interrupted if the cloud provider reclaims the machine.

To adapt this script for other CP2K simulations, replace `input_dir` with the
path to your CP2K input files and set the `sim_config_filename` accordingly.

When the simulation is complete, we terminate the machine, download the results and print a summary of the simulation as shown below.

```
Task status: Success

Timeline:
	Waiting for Input         at 17/04, 15:02:28      0.883 s
	In Queue                  at 17/04, 15:02:29      37.933 s
	Preparing to Compute      at 17/04, 15:03:07      11.183 s
	In Progress               at 17/04, 15:03:18      84.261 s
		└> 84.152 s        /opt/openmpi/5.0.6/bin/mpirun --use-hwthread-cpus cp2k.psmp H2O-64.inp
	Finalizing                at 17/04, 15:04:43      0.428 s
	Success                   at 17/04, 15:04:43      

Data:
	Size of zipped output:    86.81 KB
	Size of unzipped output:  290.38 KB
	Number of output files:   5

Estimated computation cost (US$): 0.0031 US$
```

As you can see in the "In Progress" line, the part of the timeline that represents the actual execution of the simulation, 
the core computation time of this simulation was approximately 84 seconds (around 1 min and 24 seconds).

For comparison, the same simulation takes **1 minute and 15 seconds** on a similar local machine with 16 virtual CPUs (Ryzen 7 7700X). This performance difference is expected as cloud CPUs typically have lower clock speeds compared to regular desktop processors, prioritizing energy efficiency and density over raw speed.

However, increasing the number of vCPUs on the cloud machine could improve performance.

## Scaling Up Your Simulation  
Scaling up your simulation is as simple as updating the `machine_type` parameter to a 56 vCPU machine (`c2d-highcpu-56`).

As mentioned above, running this simulation on a **16 vCPU** cloud machine was slower than on a similarly powered local computer. To improve performance, we upgraded to a **c2d-highcpu-56** instance with **56 vCPUs**, reducing the runtime to just **43 seconds** — with a slight cost increase to **$0.0058**.

| Machine Type            | vCPUs | Execution Time | Estimated Cost (USD)|
|-------------------------|-------|----------------|---------------------|
| **Local Ryzen 7 7700X** | 16    | 1 min and 15s  | N/A                 |
| **Cloud c2d-highcpu-16**| 16    | 1 min and 24s  | 0.0031              |
| **Cloud c2d-highcpu-56**| 56    | 43s            | 0.0058              | 

By leveraging the Inductiva API, you can efficiently scale your CP2K simulations
to meet your computational needs. Try different machine configurations and
optimize your workflow for faster, more cost-effective results!
