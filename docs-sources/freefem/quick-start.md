# Run Your First Simulation
This tutorial will show you how to run FreeFEM simulations using the Inductiva API. 

We will cover the `NSCaraCyl` case from the [FreeFEM Github](https://github.com/FreeFem/FreeFem-sources) to help you get started with simulations.

## Prerequisites
Download the required files [here](https://github.com/FreeFem/FreeFem-sources/blob/master/examples/mpi/NSCaraCyl.edp) and place the simulation files inside a `freefem-input-files` folder. Then, you’ll be ready to send your simulation to the Cloud.

## Running a FreeFEM Simulation
Here is the code required to run a FreeFEM simulation using the Inductiva API:

```python
"""FreeFEM example"""
import inductiva

# Instantiate machine group
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c3d-highcpu-4")

# Initialize the Simulator
freefem = inductiva.simulators.FreeFEM( \
    version="4.15")

task = freefem.run( \
	input_dir="/Path/to/freefem-input-files",
	commands=[
		"ff-mpirun -np 4 --use-hwthread-cpus NSCaraCyl.edp -cas 2 -n 40 -ndt 100 -T 5"],
	on=cloud_machine)

task.wait()
cloud_machine.terminate()

task.download_outputs()

task.print_summary()

```

In this basic example, we're using a cloud machine (`c3d-highcpu-4`) equipped with 4 virtual CPUs. 
For larger or more compute-intensive simulations, consider adjusting the `machine_type` parameter to select 
a machine with more virtual CPUs and increased memory capacity. You can explore the full range of available machines [here](https://console.inductiva.ai/machine-groups/instance-types).

> **Note**: Setting `spot=True` enables the use of [spot machines](../how-it-works/machines/spot-machines.md), which are available at substantial discounts. 
> However, your simulation may be interrupted if the cloud provider reclaims the machine.

To adapt this script for other FreeFEM simulations, replace `input_dir` with the
path to your FreeFEM input files.

When the simulation is complete, we terminate the machine, download the results and print a summary of the simulation as shown below.

```
Task status: Success

Timeline:
	Waiting for Input         at 24/09, 10:26:46      0.899 s
	In Queue                  at 24/09, 10:26:47      39.441 s
	Preparing to Compute      at 24/09, 10:27:27      5.168 s
	In Progress               at 24/09, 10:27:32      827.994 s
		└> 827.837 s       ff-mpirun -np 4 --use-hwthread-cpus NSCaraCyl.edp -cas 2 -n 40 -ndt 100 -T 5
	Finalizing                at 24/09, 10:41:20      0.559 s
	Success                   at 24/09, 10:41:20      

Data:
	Size of zipped output:    331.90 KB
	Size of unzipped output:  1.08 MB
	Number of output files:   6

Estimated computation cost (US$): 0.0097 US$
```

As you can see in the "In Progress" line, the part of the timeline that represents the actual execution of the simulation, 
the core computation time of this simulation was approximately 13 minutes and 47 seconds.

## Scaling Up Your Simulation
Running your simulation on a larger machine only requires a couple of small adjustments to your Python script.

### Changes to Apply

Update the following parameters:

* **`machine_type`**

  * Set `machine_type="c3d-highcpu-8"`
* **MPI processes (`-np`)**

  * Set `-np 8`

With these two updates, you’re doubling the computational power available to your simulation, from **4 vCPUs to 8 vCPUs**.

> **Note:**
> The final number in `machine_type` indicates the total number of vCPUs on that machine. The `-np` value specifies how many of those vCPUs your simulation will use.
>
> * `-np` can be lower than the number of vCPUs,
> * but it cannot exceed it.
>   Adjust both values to experiment with different configurations.

### Performance Comparison

The table below shows how execution time and cost change when scaling the same simulation across larger machines:

| Machine Type       | vCPUs | Execution Time | Estimated Cost (USD) |
| ------------------ | ----- | -------------- | -------------------- |
| **c3d-highcpu-4**  | 4     | 13m 47s        | 0.0097               |
| **c3d-highcpu-8**  | 8     | 8m 42s         | 0.012                |
| **c3d-highcpu-16** | 16    | 5m 31s         | 0.014                |
| **c3d-highcpu-30** | 30    | 4m 27s         | 0.021                |
| **c3d-highcpu-60** | 60    | 3m 59s         | 0.038                |

### Finding the Right Balance

With the **Inductiva API**, scaling your FreeFEM simulations is straightforward. Whether your priority is minimizing runtime or optimizing costs, experimenting with different machine configurations will help you strike the best balance for your needs.

```{banner_small}
:origin: freefem_quick_start
```