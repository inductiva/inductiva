# Run Your First Simulation
This tutorial will show you how to run REEF3D simulations using the Inductiva API. 

We will cover the `9_1 Regular Wave Propagation` use case from the FNPF Tutorials, available in the [REEF3D GitHub repository](https://github.com/REEF3D/REEF3D/tree/ed0c8d7a6110892706357f72e0404bd63034efa5), to help you get started with simulations.

## Prerequisites
Download the required files [here](https://github.com/REEF3D/REEF3D/tree/master/Tutorials/REEF3D_FNPF/9_1%20Regular%20Wave%20Propagation) and place them in a folder called `RegularWavePropagation`. Then, you’ll be ready to send your simulation to the Cloud.

## Running a REEF3D Simulation
Reef3D in the Inductiva API executes two sequential steps:
- Meshing with **DiveMESH**
- Simulation with **Reef3D**

Each step is configured with input files, `control.txt` and `ctrl.txt` respectively. Other files may be used to inform the simulator about the grid, geographical data or wave information. Reef3D has strict naming rules for each file and we recommend that users follow their [guidelines](https://reef3d.wordpress.com/user-guide/).

Here is the code required to run a REEF3D simulation using the Inductiva API:

```python
"""REEF3D example."""
import inductiva

# Allocate cloud machine on Google Cloud Platform
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c2d-highcpu-8",
    spot=True)

# Initialize the Simulator
reef3d = inductiva.simulators.REEF3D( \
    version="25.02")

# Run simulation
task = reef3d.run(input_dir="/Path/to/RegularWavePropagation",
    n_vcpu=16,
    on=cloud_machine)

# Wait for the simulation to finish and download the results
task.wait()
cloud_machine.terminate()

task.download_outputs()

task.print_summary()
```

In this basic example, we're using a cloud machine (`c2d-highcpu-8`) equipped with 8 virtual CPUs. 
For larger or more compute-intensive simulations, consider adjusting the `machine_type` parameter to select 
a machine with more virtual CPUs and increased memory capacity. You can explore the full range of available machines [here](https://console.inductiva.ai/machine-groups/instance-types).

> **Note**: Setting `spot=True` enables the use of spot machines, which are available at substantial discounts. 
> However, your simulation may be interrupted if the cloud provider reclaims the machine.

The number of virtual CPUs (`n_vcpus`) is the parameter used to configure the simulation parallelism. This value must be consistently set to the same parameter `M 10` in both the `control.txt` and `ctrl.txt` configuration files.

To adapt this script for other REEF3D simulations, replace `input_dir` with the
path to your REEF3D input files.

When the simulation is complete, we terminate the machine, download the results and print a summary of the simulation as shown below.

```
Task status: Success

Timeline:
	Waiting for Input         at 21/04, 18:27:21      1.094 s
	In Queue                  at 21/04, 18:27:22      63.809 s
	Preparing to Compute      at 21/04, 18:28:26      7.884 s
	In Progress               at 21/04, 18:28:34      20.295 s
		├> 1.061 s         /DIVEMesh/bin/DiveMESH
		└> 19.075 s        /opt/openmpi/4.1.6/bin/mpirun --use-hwthread-cpus /REEF3D/bin/REEF3D
	Finalizing                at 21/04, 18:28:54      1.283 s
	Success                   at 21/04, 18:28:56      

Data:
	Size of zipped output:    14.24 MB
	Size of unzipped output:  40.16 MB
	Number of output files:   764

Estimated computation cost (US$): 0.00051 US$
```

As you can see in the "In Progress" line, the part of the timeline that represents the actual execution of the simulation, 
the core computation time of this simulation was approximately 20 seconds.

It's that simple!
