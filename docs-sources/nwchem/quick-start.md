# Run Your First Simulation
This tutorial will show you how to run NWChem simulations using the Inductiva API. 

We will cover the `Single point Water SCF energy` use case from the examples of the [official NWChem documentation](https://nwchemgit.github.io/Sample.html) to help you get started with simulations.

## Prerequisites
1. Copy the following input file exactly as shown:
```
 start h2o 
 title "Water in 6-31g basis set" 

 geometry units au  
   O      0.00000000    0.00000000    0.00000000  
   H      0.00000000    1.43042809   -1.10715266  
   H      0.00000000   -1.43042809   -1.10715266 
 end  
 basis  
   H library 6-31g  
   O library 6-31g  
 end
 task scf
 ```
2. Create a file named `water_scf.nw` and paste the above input into it.
3. Save this file inside a folder named `SimulationFiles`.

Once these steps are complete, you’ll be ready to send your simulation to the Cloud.

## Running an NWChem Simulation
Here is the code required to run an NWChem simulation using the Inductiva API:

```python
"""NWChem example"""
import inductiva

# Allocate cloud machine on Google Cloud Platform
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c2d-highcpu-16",
    spot=True)


# Initialize the Simulator
nwchem = inductiva.simulators.NWChem(\
    version="7.2.3")

# Run simulation
task = nwchem.run(input_dir="/Path/to/SimulationFiles",
    sim_config_filename="water_scf.nw",
    n_vcpus=8,
    on=cloud_machine)

# Wait for the simulation to finish and download the results
task.wait()
cloud_machine.terminate()

task.download_outputs()

task.print_summary()
```

In this basic example, we're using a cloud machine (`c2d-highcpu-16`) equipped with 16 virtual CPUs. 
For larger or more compute-intensive simulations, consider adjusting the `machine_type` parameter to select 
a machine with more virtual CPUs and increased memory capacity. You can explore the full range of available machines [here](https://console.inductiva.ai/machine-groups/instance-types).

> **Note**: Setting `spot=True` enables the use of [spot machines](../how-it-works/machines/spot-machines.md), which are available at substantial discounts. 
> However, your simulation may be interrupted if the cloud provider reclaims the machine.

To adapt this script for other NWChem simulations, replace `input_dir` with the
path to your NWChem input files and set the `sim_config_filename` accordingly.

When the simulation is complete, we terminate the machine, download the results and print a summary of the simulation as shown below.

```
Task status: Success

Timeline:
	Waiting for Input         at 21/04, 15:23:46      0.815 s
	In Queue                  at 21/04, 15:23:47      36.562 s
	Preparing to Compute      at 21/04, 15:24:23      3.109 s
	In Progress               at 21/04, 15:24:27      2.355 s
		└> 2.253 s         /opt/openmpi/4.1.6/bin/mpirun --use-hwthread-cpus --np 8 nwchem water_scf.nw
	Finalizing                at 21/04, 15:24:29      0.39 s
	Success                   at 21/04, 15:24:29      

Data:
	Size of zipped output:    19.64 KB
	Size of unzipped output:  824.11 KB
	Number of output files:   10

Estimated computation cost (US$): 0.00020 US$
```

As you can see in the "In Progress" line, the part of the timeline that represents the actual execution of 
the simulation, the core computation time of this simulation was approximately 2 seconds.

```{banner_small}
:origin: nwchem
```

