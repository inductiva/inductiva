# Run Your First Simulation
This tutorial will show you how to run SFINCS simulations using the Inductiva API. 

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

## Running an SFINCS Simulation
Here is the code required to run a SCHISM simulation using the Inductiva API:

```python
"""SFINCS example"""
import inductiva

# Allocate cloud machine on Google Cloud Platform
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c2d-highcpu-8",
	spot=True)

# Initialize the Simulator
sfincs = inductiva.simulators.SFINCS(\
    version="2.2.1")

# Run simulation
task = sfincs.run(\
    input_dir="/Path/to/Models/HighReliefArea/SFINCS/InputFiles",
    on=cloud_machine)

# Wait for the simulation to finish and download the results
task.wait()
cloud_machine.terminate()

task.download_outputs()

task.print_summary()
```

> **Note**: `spot` machines are a lot cheaper but may be terminated by the provider if necessary.

To adapt this script for other SFINCS simulations, replace `input_dir` with the
path to your SFINCS input files.

When the simulation is complete, we terminate the machine, download the results and print a summary of the simulation as shown below.

```
Task status: Success

Timeline:
	Waiting for Input         at 04/06, 14:53:34      3.135 s
	In Queue                  at 04/06, 14:53:37      38.267 s
	Preparing to Compute      at 04/06, 14:54:16      2.278 s
	In Progress               at 04/06, 14:54:18      1823.975 s
		└> 1823.822 s      sfincs
	Finalizing                at 04/06, 15:24:42      2.335 s
	Success                   at 04/06, 15:24:44      

Data:
	Size of zipped output:    74.44 MB
	Size of unzipped output:  88.29 MB
	Number of output files:   7

Estimated computation cost (US$): 0.024 US$
```

As you can see in the “In Progress” line, the part of the timeline that represents the actual execution of the simulation, the core computation time of this simulation was approximately 1824 seconds (30 minutes and 24 seconds).

<br>

> 💻 Want to run SFINCS directly from Deltares' repository? Check out our [tutorial](https://inductiva.ai/guides/how-it-works/bring-your-own-software/run-sfincs-directly-from-deltares-repository) to get started.