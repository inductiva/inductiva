# Run a Single Simulation
First, we will run a single OpenFAST simulation using the `5MW_OC4Semi_WSt_WavesWN` case. This should be straightforward as all 
the necessary input files are already prepared.

## Code Overview
The code required to run an OpenFAST simulation using the Inductiva API is always the same, we just need to adapt it for this use case, as shown below.

```python
import inductiva

# Allocate Cloud machine on Google Cloud Platform
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

# Wait for the simulation to finish and download the results
task.wait()
cloud_machine.terminate()

task.download_outputs()

task.print_summary()
```

When the simulation is complete, we terminate the machine, download the results and print a summary of the simulation as follows.

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

## Performance and Cost Analysis
Given that OpenFAST does not benefit from multiple CPU cores, we chose the `n2-highcpu-2` virtual machine (VM) with 2 virtual CPUs (equivalent to 1 physical core). 
This is one of the cheapest options on Google Cloud, costing just US$0.0081 per hour in spot mode.

To demonstrate that OpenFAST does not scale with the number of cores, we also ran the same simulation on a number of better machines. For a detailed breakdown, check out our [Benchmarks](https://inductiva.ai/guides/openfast/Benchmarks) section.

Here are the results:
| Machine       | Number of VCPUs | Execution time | Cost |
|---------------|-----------------|----------------|------|
| n2-highcpu-2  | 2               |34.1 seconds  |0.00013 US$|
| n2-highcpu-4  | 4               |35.0 seconds  |0.00031 US$|
| n2-highcpu-8  | 8               |30.4 seconds  |0.00034 US$|
| n2-highcpu-16 | 16              |30.9 seconds    |0.00065 US$|
| n2-highcpu-32 | 32              |30.2 seconds  |0.0012 US$|

The execution time remains almost the same on all machines, regardless of the number of virtual CPUs.

Additionally, running the simulation on the `n2-highcpu-2` VM proves to be extremely cost-efficient, with a total cost of only 0.00013 US$.

In the next part of this tutorial, we'll take things to the next level by running dozens of OpenFAST simulations in parallel on Inductiva, demonstrating 
the true power of cloud-based scalability. Stay tuned!