# Run Your First Simulation
This tutorial will show you how to run OpenFOAM simulations using the Inductiva API. 

We will cover the `motorBike` use case from the [OpenFOAM Foundation GitHub repository](https://github.com/OpenFOAM/OpenFOAM-8/tree/version-8/tutorials), to help you get started with simulations.

## Prerequisites
Download the required files [here](https://github.com/OpenFOAM/OpenFOAM-8/tree/version-8/tutorials/incompressible/simpleFoam/motorBike) and place them in a folder called `SimulationFiles`.

Now, open the file `system/controlDict` and change the following lines:

```diff
-numberOfSubdomains 6;
+numberOfSubdomains 8;

method          hierarchical;
// method          ptscotch;

simpleCoeffs
{
    n               (4 1 1);
    delta           0.001;
}

hierarchicalCoeffs
{
-    n               (3 2 1);
+    n               (4 2 1);
    delta           0.001;
    order           xyz;
}
```

This change sets the simulation to run with 8 processes, matching the number of
**physical cores** available on a `c2d-highcpu-16` machine (that will be used for this tutorial).

> **Note**: In OpenFOAM-Foundation v8, the `runParallel` function limits the simulation to the number of **physical** cores on the machine. That’s why we configure 8 subdomains: to fully utilize all physical cores available on this machine.
> Learn more about this [here](faq.md#6-why-does-my-simulation-keep-failing-with-there-are-not-enough-slots-available-even-though-my-machine-has-enough-resources).

## Running an OpenFOAM Simulation
Here is the code required to run an OpenFOAM simulation using the Inductiva API:

```python
"""OpenFOAM example"""
import inductiva

# Allocate cloud machine on Google Cloud Platform
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c2d-highcpu-16",
	# 1 thread per physical core
	threads_per_core=1,
	spot=True)

# Initialize the Simulator
OpenFOAM = inductiva.simulators.OpenFOAM( \
    version="8",
	distribution="foundation")

# Run simulation
task = OpenFOAM.run(input_dir="/Path/to/SimulationFiles",
    shell_script="./Allrun",
	# add simulation to a project
	project="openfoam-quick-start",
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

> **Note**: Setting `spot=True` enables the use of spot machines, which are available at substantial discounts. 
> However, your simulation may be interrupted if the cloud provider reclaims the machine.

To adapt this script for other OpenFOAM simulations, replace `input_dir` with the
path to your OpenFOAM input files and set the OpenFOAM `distribution` and `version` accordingly.

We run the simulation using the `run` method, specifying the `shell_script` that handles the execution process.

When the simulation is complete, we terminate the machine, download the results and print a summary of the simulation as shown below.

```
inductiva tasks info qc1b1yctqhxnf0knjpdswxhyw

Task status: Success

Timeline:
	Waiting for Input         at 14/07, 10:51:24      0.796 s
	In Queue                  at 14/07, 10:51:25      38.967 s
	Preparing to Compute      at 14/07, 10:52:04      4.014 s
	In Progress               at 14/07, 10:52:08      112.572 s
		└> 112.388 s       bash ./Allrun
	Finalizing                at 14/07, 10:54:00      1.752 s
	Success                   at 14/07, 10:54:02      

Data:
	Size of zipped output:    253.07 MB
	Size of unzipped output:  346.32 MB
	Number of output files:   600

Estimated computation cost (US$): 0.0030 US$

Go to https://console.inductiva.ai/tasks/qc1b1yctqhxnf0knjpdswxhyw for more details.
```

As you can see in the "In Progress" line, the part of the timeline that represents the actual execution of the simulation, 
the core computation time of this simulation was approximately 1 minutes and 52 seconds.

Here’s a clearer and more polished version of your section:

---

## Scaling Up the Simulation

After successfully running the simulation on 8 physical cores, let’s scale it up
to 16 physical cores. This corresponds to a machine with 32 virtual CPUs.

### Update the Machine Configuration

Modify your script to use a more powerful machine by setting `machine_type` to
one with 32 vCPUs (i.e., 16 physical cores).

```python
cloud_machine = inductiva.resources.MachineGroup(
    provider="GCP",
    machine_type="c2d-highcpu-32",
    # 1 thread per physical core
	threads_per_core=1,
    spot=True
)
```

### Adjust the Simulation Configuration

Update the `system/controlDict` file to divide the simulation into 16 sub-domains:

```diff
-numberOfSubdomains 8;
+numberOfSubdomains 16;

method          hierarchical;
// method          ptscotch;

simpleCoeffs
{
    n               (4 1 1);
    delta           0.001;
}

hierarchicalCoeffs
{
-    n               (3 2 1);
+    n               (8 2 1);
    delta           0.001;
    order           xyz;
}
```

Below are the results from running the simulation on both 8 and 16 physical cores:

| Machine Type       | Sub-domains | Duration (hh:mm:ss) | Cost (USD) |
|--------------------|-------------|----------------------|------------|
| c2d-highcpu-16     | 8           | 1 min 52 sec        | 0.0030 US$  |
| c2d-highcpu-32     | 16          | 1 min 8 sec         | 0.0037 US$  |

As shown by the results, we managed to run the simulation **1.65 times** faster
with only an increase in cost of **1.23 times**.


<div class="cta-bar">
  <div class="cta-text">
    <strong>Start running simulations seamlessly!</strong> You have $5 in <strong>free</strong> credits, no credit card required.
  </div>
  <button  onclick="window.open('https://console.inductiva.ai/', '_blank')" target="_blank" class="cta-button">Sign In</button>
</div>

