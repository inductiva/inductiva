# Run your UnSWAN (UNstructured mesh SWAN) simulation
This tutorial will show you how to run UnSWAN simulations using the Inductiva API. 

We will cover the `Haringvliet field case` from the [official SWAN documentation](https://swanmodel.sourceforge.io/download/download.htm).

## Prerequisites
Download the required files [here](https://swanmodel.sourceforge.io/download/zip/f32harin.tar.gz) and place them in a folder called `f32harin`. Then, you’ll be ready to send your simulation to the Cloud.

## Running an UnSWAN Simulation
Here is the code required to run a UnSWAN simulation using the Inductiva API:

```python
"""UnSWAN example"""
import inductiva

# Allocate cloud machine on Google Cloud Platform
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c2d-highcpu-4",
	spot=True)

# Initialize the Simulator
swan = inductiva.simulators.SWAN(\
    version="41.51")

# Run simulation
task = swan.run(input_dir="/Path/to/f32harin",
    sim_config_filename="f32har01.swn",
	# run unswan
	command="unswan",
    on=cloud_machine)

# Wait for the simulation to finish and download the results
task.wait()
cloud_machine.terminate()

task.download_outputs()

task.print_summary()
```

In this basic example, we're using a cloud machine (`c2d-highcpu-4`) equipped with 4 virtual CPUs. 
For larger or more compute-intensive simulations, consider adjusting the `machine_type` parameter to select 
a machine with more virtual CPUs and increased memory capacity. You can explore the full range of available machines [here](https://console.inductiva.ai/machine-groups/instance-types).

> **Note**: Setting `spot=True` enables the use of [spot machines](../how-it-works/machines/spot-machines.md), which are available at substantial discounts. 
> However, your simulation may be interrupted if the cloud provider reclaims the machine.

To adapt this script for other SWAN simulations, replace `input_dir` with the
path to your SWAN input files and set the `sim_config_filename` accordingly.

When the simulation is complete, we terminate the machine, download the results and print a summary of the simulation as shown below.

```
Task status: Success

Timeline:
	Waiting for Input         at 17/09, 14:36:59      0.945 s
	In Queue                  at 17/09, 14:37:00      43.861 s
	Preparing to Compute      at 17/09, 14:37:44      4.474 s
	In Progress               at 17/09, 14:37:48      28.266 s
		└> 28.093 s        unswanrun -input f32har01.swn -omp 4
	Finalizing                at 17/09, 14:38:16      0.504 s
	Success                   at 17/09, 14:38:17      

Data:
	Size of zipped output:    345.55 KB
	Size of unzipped output:  536.86 KB
	Number of output files:   10

Estimated computation cost (US$): 0.00024 US$
```

As you can see in the "In Progress" line, the part of the timeline that represents the actual execution of the simulation, 
the core computation time of this simulation was 28.3 seconds.

```{banner_small}
:origin: swan-unswan-tutorial
```
