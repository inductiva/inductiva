# Run a Single Simulation
Let’s begin by running a single OpenFOAM simulation using the `occDrivAerStaticMesh` case.
All the required input files are already prepared, as outlined in the previous section.

## Code Overview
The Python code required to run a OpenFOAM simulation using the Inductiva API follows a consistent structure. We adapted it for this specific use case, as shown below.

```python
import inductiva

# Allocate cloud machine on Google Cloud Platform
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c3d-highcpu-180",
    data_disk_gb=100,
    spot=True)

# Initialize OpenFOAM stack
openfoam = inductiva.simulators.OpenFOAM(
    version="2412",
    distribution="esi"
    )

task = openfoam.run( \
    input_dir="/Path/to/openfoam-occDrivAerStaticMesh",
    shell_script="./Allrun",
    on=cloud_machine)

# Wait for the simulation to finish and download the results
task.wait()
cloud_machine.terminate()

task.download_outputs()

task.print_summary()
```

When the simulation is complete, we terminate the machine, download the results and print a summary of the simulation 
as shown.

```
inductiva tasks info 6j356wymfn7t614fncoff2dnf
Task status: Success

Timeline:
	Waiting for Input         at 09/07, 09:25:15      1933.174 s
	In Queue                  at 09/07, 09:57:28      4.302 s
	Preparing to Compute      at 09/07, 09:57:32      53.658 s
	In Progress               at 09/07, 09:58:26      51074.226 s
		└> 51074.226 s        bash Allrun
	Finalizing                at 10/07, 00:09:40      23.903 s
	Success                   at 10/07, 00:10:04      

Data:
	Size of zipped output:    5.51 GB
	Size of unzipped output:  10.72 GB
	Number of output files:   3276

Estimated computation cost (US$): 24.31 US$

Go to https://console.inductiva.ai/tasks/6j356wymfn7t614fncoff2dnf for more details.
```

As you can see in the “In Progress” line, the part of the timeline that represents the actual execution of the simulation, 
the core computation time of this simulation was approximately 14 hours and 13 minutes.

Next, we’ll show you how to scale this same simulation across multiple machines and significantly 
reduce simulation time by a factor of **1.8 times**!

Stay tuned!