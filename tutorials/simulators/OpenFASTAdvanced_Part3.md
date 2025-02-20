---
orphan: true
---

# Running our Simulation

Now that you have all the necessary files set up from the previous section,
you're ready to run your first simulation. In this part of the tutorial, we will
execute a single OpenFAST simulation using the `5MW_OC4Semi_WSt_WavesWN` case.
The process should be straightforward, as all required input files are already
in place.

We'll go step by step to ensure everything runs smoothly, so let's get started!


## Overview

For you to run this simulation you can execute the following python code.

```python
import inductiva

# Allocate cloud machine
cloud_machine = inductiva.resources.MachineGroup(
    provider="GCP",
    machine_type="n2-highcpu-2",
    spot=True)

# Initialize OpenFAST simulator
openfast = inductiva.simulators.OpenFAST(
    version="4.0.2")

# Run simulation
task = openfast.run(
    input_dir="input_files",
    commands=[
        "openfast 5MW_OC4Semi_WSt_WavesWN/"
        "5MW_OC4Semi_WSt_WavesWN.fst"],
    on=cloud_machine)

task.wait()
cloud_machine.terminate()
task.download_outputs()
task.print_summary()

```

### Simulation Results

After running the simulation, here is the task summary:

```
inductiva tasks info hk8fmb4dif05c5vf5saurocxp

Task status: Success

Timeline:
	Waiting for Input         at 20/02, 14:24:28      6.487 s
	In Queue                  at 20/02, 14:24:35      19.415 s
	Preparing to Compute      at 20/02, 14:24:54      1.842 s
	In Progress               at 20/02, 14:24:56      34.265 s
		└> 34.123 s        openfast 5MW_OC4Semi_WSt_WavesWN/5MW_OC4Semi_WSt_WavesWN.fst
	Finalizing                at 20/02, 14:25:30      1.085 s
	Success                   at 20/02, 14:25:31      

Data:
	Size of zipped output:    16.20 MB
	Size of unzipped output:  32.71 MB
	Number of output files:   83

Estimated computation cost (US$): 0.00013 US$
```

### Performance and Cost Analysis

Since OpenFAST does not benefit from multiple CPU cores, we selected the
`n2-highcpu-2` virtual machine (VM) with 2 virtual CPUs (equivalent to 1 physical
core). This is one of the cheapest options on Google Cloud, costing just
0.0081 US$/hour in spot mode.

To verify that OpenFAST does not scale with the number of cores, we also ran the
same simulation on a n2-highcpu-8 VM with 8 virtual CPUs:

```
inductiva tasks info s3mph8i4cbq0ute0011nqisfe

Task status: Success

Timeline:
	Waiting for Input         at 20/02, 14:25:15      7.525 s
	In Queue                  at 20/02, 14:25:23      4.923 s
	Preparing to Compute      at 20/02, 14:25:28      1.831 s
	In Progress               at 20/02, 14:25:29      30.409 s
		└> 30.274 s        openfast 5MW_OC4Semi_WSt_WavesWN/5MW_OC4Semi_WSt_WavesWN.fst
	Finalizing                at 20/02, 14:26:00      1.12 s
	Success                   at 20/02, 14:26:01      

Data:
	Size of zipped output:    16.20 MB
	Size of unzipped output:  32.71 MB
	Number of output files:   83

Estimated computation cost (US$): 0.00035 US$
```

The results confirm what we knew to be true: the execution time (In Progress
duration) remains nearly the same on a machine with four times the number of
virtual CPUs.

Furthermore, running the simulation on the n2-highcpu-2 VM proves to be
extremely cost-efficient, with a total cost of just 0.00013 US$.

In the next part of this tutorial, we'll take things to the next level by
running dozens of OpenFAST simulations in parallel using Inductiva,
demonstrating the true power of cloud-based scalability. Stay tuned!

[Running 40 Simulations - Templating](OpenFASTAdvanced_Part4.md)