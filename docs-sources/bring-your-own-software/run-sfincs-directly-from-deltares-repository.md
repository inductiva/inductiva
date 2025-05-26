# Run SFINCS Directly from Deltares’ Repository
The [SFINCS](https://www.deltares.nl/en/software-and-data/products/sfincs) (Super-Fast INundation of CoastS) model is a reduced-complexity simulation 
engine developed by [Deltares](https://www.deltares.nl/en). It is designed 
to simulate compound coastal flooding with high computational efficiency 
while maintaining reliable accuracy.

This tutorial shows you how to run SFINCS through the Inductiva API using 
Deltares' publicly available [container image](https://hub.docker.com/r/deltares/sfincs-cpu) directly from their repository.

The use case covered in this tutorial is featured in the following publication: *[Camila Gaido Lasserre, Kees Nederhoff, Curt D. Storlazzi, Borja G. Reguero, Michael W. Beck (2024). "Improved efficient physics-based computational modeling of regional wave-driven coastal flooding for reef-lined coastlines." Ocean Modelling.](https://www.sciencedirect.com/science/article/pii/S1463500324000453#refdata001)*

## Prerequisites
Download the use case [here](https://zenodo.org/records/10805615) and copy 
the `sfincs_netbndbzsbzifile.nc` file from `Models/HighReliefArea/SFINCS/BoundaryConditions/rp_000_SLR000` to the `Models/HighReliefArea/SFINCS/InputFiles` folder.

The `InputFiles` folder should contain the following files:

```
sfincs.dep			sfincs.man
sfincs.gms			sfincs.msk
sfincs.ind			sfincs.obs
sfincs.inp			sfincs.obs.bak
sfincs.inp.bak		sfincs_netbndbzsbzifile.nc
```

## Run Simulations using SFINCS Docker Container
You can reference the Deltares’ Docker image directly in your simulation script using the `inductiva.simulators.CustomImage`class, as shown below:

```python
import inductiva

cloud_machine = inductiva.resources.MachineGroup(
	"c2d-highcpu-8",
	data_disk_gb=100)

custom_simulator = inductiva.simulators.CustomImage(
	container_image="deltares/sfincs-cpu:latest")

input_dir = "/Path/to/Models/HighReliefArea/SFINCS/InputFiles"

task = custom_simulator.run(
	input_dir=input_dir,
	commands=["sfincs"],
	on=cloud_machine)
```

The `CustomImage` simulator enables you to run simulations using any Docker image of your choice by specifying the `container_image` parameter. 

When the simulation is complete, we terminate the machine, download the results and print a summary of the simulation as shown below.

```
Task status: Success

Timeline:
	Waiting for Input         at 07/05, 15:59:41      2.927 s
	In Queue                  at 07/05, 15:59:43      34.238 s
	Preparing to Compute      at 07/05, 16:00:18      22.302 s
	In Progress               at 07/05, 16:00:40      1443.74 s
		└> 1443.575 s      sfincs
	Finalizing                at 07/05, 16:24:44      2.737 s
	Success                   at 07/05, 16:24:47      

Data:
	Size of zipped output:    74.44 MB
	Size of unzipped output:  88.29 MB
	Number of output files:   7

Estimated computation cost (US$): 0.028 US$
```

As you can see in the “In Progress” line, the part of the timeline that represents the actual execution of the simulation, the core computation time of this simulation was approximately 1443.7 seconds (24 minutes and 44 seconds).

## Scaling Up Your Simulation
Scaling up your simulation is as simple as updating the `machine_type` parameter to a 16 vCPU machine (c2d-highcpu-16).

By increasing the number of virtual CPUs, we’ve reduced the processing time from 25 minutes and 31 seconds to 18 minutes and 21 seconds.

Here are the results of running the same simulation on a few machines:

|   Machine Type  | Virtual CPUs |     Execution Time     |   Estimated Cost   |
|:---------------:|:------------:|:----------------------:|:--------:|
|  c2d-highcpu-8  |       8      | 24 minutes and 44 seconds | 0.028 US$ |
|  c2d-highcpu-16 |      16      |  18 minutes and 21 seconds | 0.037 US$ |
|  c2d-highcpu-32 |      32      |  12 minutes and 57 seconds  | 0.048 US$ |
|  c2d-highcpu-56 |      56     |   9 minutes and 47 seconds  | 0.064 US$ |
|  c3d-highcpu-8 |      8    |    24 minutes and 38 seconds  | 0.039 US$ |
|  c3d-highcpu-16 |      16    |  19 minutes and 23 seconds  |  0.057 US$ |
|  c3d-highcpu-30 |      30   |  12 minutes and 55 seconds  | 0.065 US$ |
|  c3d-highcpu-60 |      60  |   9 minutes and 45 seconds  | 0.091 US$ |

It's that simple!