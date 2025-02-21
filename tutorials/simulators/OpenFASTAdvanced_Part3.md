---
orphan: true
---

# Running the Simulation using Inductiva API

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
# Versions available: 3.5.2 and 4.0.2 
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
same simulation on a set of better machines and here are the results:

| Machine       | Number of VCPUs | Execution time | Cost |
|---------------|-----------------|----------------|------|
| n2-highcpu-2  | 2               |34.1 seconds  |0.00013 US$|
| n2-highcpu-4  | 4               |35.0 seconds  |0.00031 US$|
| n2-highcpu-8  | 8               |30.4 seconds  |0.00034 US$|
| n2-highcpu-16 | 16              |30.9 seconds    |0.00065 US$|
| n2-highcpu-32 | 32              |30.2 seconds  |0.0012 US$|

The results confirm what we knew to be true: the execution time remains nearly
the same on all machines regardless of the number of virtual CPUs.

Furthermore, running the simulation on the n2-highcpu-2 VM proves to be
extremely cost-efficient, with a total cost of just 0.00013 US$.

In the next part of this tutorial, we'll take things to the next level by
running dozens of OpenFAST simulations in parallel using Inductiva,
demonstrating the true power of cloud-based scalability. Stay tuned!

[Running 50 Simulations in parallel - Templating](OpenFASTAdvanced_Part4.md)