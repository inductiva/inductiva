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
2. Create a file named `Water_SCF.nw` and paste the above input into it.
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
    machine_type="c3d-standard-16",
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

> **Note**: `spot` machines are a lot cheaper but may be terminated by the provider if necessary.

To adapt this script for other NWChem simulations, replace `input_dir` with the
path to your NWChem input files and set the `sim_config_filename` accordingly.

When the simulation is complete, we terminate the machine, download the results and print a summary of the simulation as shown below.

```
Task status: Success

Timeline:
	Waiting for Input         at 10/04, 15:46:51      0.768 s
	In Queue                  at 10/04, 15:46:51      37.174 s
	Preparing to Compute      at 10/04, 15:47:29      2.224 s
	In Progress               at 10/04, 15:47:31      3.16 s
		└> 3.06 s          /opt/openmpi/4.1.6/bin/mpirun --use-hwthread-cpus --np 8 nwchem water_scf.nw
	Finalizing                at 10/04, 15:47:34      0.483 s
	Success                   at 10/04, 15:47:34      

Data:
	Size of zipped output:    N/A until task ends
	Size of unzipped output:  824.11 KB
	Number of output files:   11

Estimated computation cost (US$): 0.00031 US$
```

As you can see in the "In Progress" line, the part of the timeline that represents the actual execution of the simulation, 
the core computation time of this simulation was approximately 3 seconds.

It's that simple!
