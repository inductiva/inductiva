# Run Your First Simulation
This tutorial will show you how to run Elmer simulations using the Inductiva API. 

We will cover the `Transcranial Magnetic Stimulation (TMS) benchmark` from the [Elmer-LinSys Github](https://github.com/ElmerCSC/elmer-linsys) to help you get started with simulations.

## Prerequisites
Download the required files [here](https://github.com/ElmerCSC/elmer-linsys/tree/main/Magnetostatics/TMS) and place the simulation files inside a `TMS` folder. Then, you’ll be ready to send your simulation to the Cloud.

## Running a Elmer Simulation
Here is the code required to run a Elmer simulation using the Inductiva API:

```python
"""Elmer example"""
import inductiva

# Instantiate machine group
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c3d-highcpu-180")

# Initialize the Simulator
elmer = inductiva.simulators.Elmer( \
    version="9.0")

# Run simulation with config files in the input directory
task = elmer.run( \
    input_dir="/Path/to/TMS",
    commands=[
		#Generating the mesh
        "python3 gen_mesh.py",
        #Running ElmerGrid on the mesh.
        #Tells ElmerGrid to convert a Gmsh mesh (.msh)
		#into Elmer format (ElmerSolver mesh files).
        "ElmerGrid 14 2 figure_8.msh -autoclean -metisbc -halo",
        #Converts the Elmer mesh files to run with
		#180 vcpus
        f"ElmerGrid 2 2 figure_8 -partdual -metiskway 180",
		#Runs ElmerSolver mpi version
        "ElmerSolver_mpi case.sif"
    ],
    on=cloud_machine)

task.wait()
cloud_machine.terminate()

task.download_outputs()

task.print_summary()

```

In this basic example, we're using a cloud machine (`c3d-highcpu-180`) equipped with 180 virtual CPUs. 
For larger or more compute-intensive simulations, consider adjusting the `machine_type` parameter to select 
a machine with more virtual CPUs and increased memory capacity. You can explore the full range of available machines [here](https://console.inductiva.ai/machine-groups/instance-types).

> **Note**: Setting `spot=True` enables the use of [spot machines](../how-it-works/machines/spot-machines.md), which are available at substantial discounts. 
> However, your simulation may be interrupted if the cloud provider reclaims the machine.

To adapt this script for other Elmer simulations, replace `input_dir` with the
path to your Elmer input files.

When the simulation is complete, we terminate the machine, download the results and print a summary of the simulation as shown below.

```
Task status: Success

Timeline:
	Waiting for Input         at 19/09, 09:51:26      0.732 s
	In Queue                  at 19/09, 09:51:26      45.916 s
	Preparing to Compute      at 19/09, 09:52:12      4.456 s
	In Progress               at 19/09, 09:52:17      1055.303 s
		├> 48.144 s        python3 gen_mesh.py
		├> 2.077 s         ElmerGrid 14 2 figure_8.msh -autoclean -metisbc -halo
		├> 8.096 s         ElmerGrid 2 2 figure_8 -partdual -metiskway 180
		└> 996.628 s       mpirun -np 180 --use-hwthread-cpus ElmerSolver_mpi case.sif
	Finalizing                at 19/09, 10:09:52      7.743 s
	Success                   at 19/09, 10:10:00      

Data:
	Size of zipped output:    1.40 GB
	Size of unzipped output:  2.54 GB
	Number of output files:   1165

Estimated computation cost (US$): 0.50 US$
```

As you can see in the "In Progress" line, the part of the timeline that represents the actual execution of the simulation, 
the core computation time of this simulation was approximately 17 minutes and 35 seconds.

```{banner_small}
:origin: elmer_quick_start
```