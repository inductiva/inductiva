# 3-D dam break impacting an elastic plate

To evaluate the performance of DualSPHysics across different GPUs, we will run a benchmark simulation based on one presented in the paper
[A fluid–structure interaction model for free-surface flows and flexible structures using smoothed particle hydrodynamics on a GPU](https://www.sciencedirect.com/science/article/pii/S0889974621000955?via%3Dihub).
Among the various simulations discussed, we focus on the `3-D dam break impacting an elastic plate` scenario.

## Simulation Overview

This simulation showcases the well known dam break problem, where a column of
water is released from a height. The catch here is the fact that there is an elastic
plate that is placed on the other side of the simulation domain.

More about this simulation can be found in section `4.5. 3-D dam break impacting an elastic plate`
of the referenced paper.

<p align="center"><img src="./_static/dam_break_elastic.gif" alt="Visualization created with Blender." width="700"></p>

## Prerequisites

### Download the simulation files

### Download the DualSPHysics package
Download the required files from `DualSPHysics_v5.4.2.zip` [here](https://dual.sphysics.org/downloads/).
You will be working within this directory and writing the Inductiva Python script there.

### Update the simulation script of the `examples/flexstruc/01_DamBreak` case

Before running the simulation, you will need to adjust the simulation script located
at `examples/flexstruc/01_DamBreak/xCaseDambreak3D_FSI_linux64_GPU.sh`.

Make the following adjustments:
1. **Update the `dirbin` variable:**
   Modify the `xCaseDambreak3D_FSI_linux64_GPU.sh` script to point to the correct binaries directory:
   ```bash
   export dirbin=/DualSPHysics/bin/linux/
   ```
2. **Remove user input prompt:**
   To enable automated execution, delete the final line in the script, which waits for user input:
   ```bash
   read -n1 -r -p "Press any key to continue..." key
   ```

These modifications will prepare the script for seamless automated execution.


## Running Your Simulation
Here is the code required to run a DualSPHysics simulation using the Inductiva API:

```python
"""DualSPHysics example."""
import inductiva

# Allocate cloud machine on Google Cloud Platform
cloud_machine = inductiva.resources.MachineGroup(
    provider="GCP",
    machine_type="g2-standard-16",
    spot=True,
    data_disk_gb=200)

# Initialize the Simulator
dualsphysics = inductiva.simulators.DualSPHysics( \
    version="5.4.1")

# Run simulation
task = dualsphysics.run(
    input_dir="Path/to/01_DamBreak",
    shell_script="xCaseDambreak3D_FSI_linux64_GPU.sh",
    on=cloud_machine)

# Wait for the simulation to finish and download the results
task.wait()
task.download_outputs()

cloud_machine.terminate()

task.print_summary()
```

This simulation runs on a `g2-standard-16` machine on spot mode, which has 16 virtual CPUs,
1 nvidia-l4 GPU and a 200 GB data disk.

> **Note**: `spot` machines are a lot cheaper but may be terminated by the provider if necessary.

When the simulation is complete, we terminate the machine, download the results and print a summary of the simulation as shown below.

```
Task status: Success

Timeline:
	Waiting for Input         at 12/05, 15:47:41      0.836 s
	In Queue                  at 12/05, 15:47:42      58.49 s
	Preparing to Compute      at 12/05, 15:48:40      13.809 s
	In Progress               at 12/05, 15:48:54      5750.479 s
		└> 5750.347 s      bash xCaseDambreak3D_FSI_linux64_GPU.sh
	Finalizing                at 12/05, 17:24:44      90.906 s
	Success                   at 12/05, 17:26:15      

Data:
	Size of zipped output:    15.06 GB
	Size of unzipped output:  31.94 GB
	Number of output files:   1029

Estimated computation cost (US$): 0.67 US$

Go to https://console.inductiva.ai/tasks/mmxoeshep1t2j2rcootdrrx5s for more details.
```

As you can see in the "In Progress" line, the part of the timeline that 
represents the actual execution of the simulation, the core computation time 
of this simulation was approximately 1 hour and 36 minutes (5750 seconds).

## Testing different GPUS

As it was stated before, the simulation ran on an Nvidia L4 GPU. What if we want
to speed up the simulation by using a more powerful GPU? We can do that by
changing the `machine_type` parameter in the code above to `a2-highgpu-1g` or
`a3-highgpu-1g`. This machines have 1 Nvidia A100 and 1 Nvidia H100 GPU respectively.

**NOTE:** you also have to add `zone="europe-west4-b"` to the machine creation due to
the fact that the A100 and H100 GPUs are not yet available on our default zone.


Here are the results of running the same simulation on those machines:

|  Machine Type  | GPU         |Execution Time          | Estimated Cost |
|:--------------:|:-----------:|:----------------------:|:--------------:|
|  g2-standard-16| Nvidia L4   | 1 hour 36 minutes      | 0.66 US$    |
|  a2-highgpu-1g | Nvidia A100 | 1 hour 2 minutes       | 0.61 US$    |
|  a3-highgpu-1g | Nvidia H100 | 34 minutes 33 seconds  | 1.55 US$    |


Keep reading to learn how to [visualize the results of your simulation using
Paraview or Blender](./visualization/index).
