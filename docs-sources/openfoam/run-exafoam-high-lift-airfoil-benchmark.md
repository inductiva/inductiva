# Running the 30P30N Micro-Benchmark (MB9) from ExaFOAM

In this tutorial, we’ll walk you through how to run the **30P30N micro-benchmark**
from the ExaFOAM suite using the **Inductiva API**. This is an advanced OpenFOAM
case that requires substantial computational power, making it a perfect showcase
for cloud-based simulation.

## 🧭 Objective

The goal of this tutorial is to demonstrate how to run the `30P30N` benchmark
from the [official ExaFOAM documentation](https://exafoam.eu/benchmarks/) using
cloud computing resources provided by Inductiva.

## 📦 Prerequisites

Before running the simulation, you need to download the benchmark files and make
a few small adjustments.

### 1. Download the Benchmark Case

Download the case files from the ExaFOAM repository:
[Download files](https://develop.openfoam.com/committees/hpc/-/tree/6b9c8438ec54961d2968977f84a47a5f9c5be4b6/compressible/rhoPimpleFoam/LES/highLiftConfiguration)

Place them in a folder named `highLiftConfiguration`.

Your folder structure should look like this:

```bash
ls -lasgo highLiftConfiguration
total 128
 0 drwxrwxr-x@ 13     416 Apr 10 11:13 .
 0 drwx------@  5     160 Apr 23 08:35 ..
24 -rw-r--r--@  1   10244 Jun 12 15:58 .DS_Store
 0 drwxrwxr-x@ 11     352 Apr  8 11:27 0.orig
 8 -rwxr-xr-x@  1     626 Apr  8 11:27 Allclean
16 -rwxr-xr-x@  1    7017 Jun 12 15:41 Allrun
 8 -rw-rw-r--@  1     991 Apr  8 11:27 COPYING
48 -rw-rw-r--@  1   21547 Apr  8 11:27 README.md
 0 -rw-rw-r--@  1       0 Apr  8 11:27 case.foam
 0 drwxrwxr-x@  5     160 Apr  8 11:27 constant
 0 drwxrwxr-x@ 13     416 Apr  8 11:27 figures
 0 drwxrwxr-x@ 28     896 Apr 10 11:13 system
24 -rw-rw-r--@  1   11399 Apr  8 11:27 thumbnail.png
```

### 2. Apply Modifications

* Open the `Allrun` script and **replace**:

  ```bash
  parEx="mpirun -np $nProcs"
  ```

  with:

  ```bash
  parEx="mpirun -use-hwthread-cpus -np $nProcs"
  ```

* Edit `highLiftConfiguration/system/include/caseDefinition` and set:

  ```bash
  nCores 360;
  ```

> **Note**: The `-use-hwthread-cpus` flag enables all available virtual CPUs on the machine for optimal performance.

## 🚀 Running the Simulation

With everything set up, you can now use the Inductiva API to run the simulation on a high-performance cloud machine.

Here’s the full Python script:

```python
"""Run ExaFOAM MB9 benchmark with Inductiva API"""
import inductiva

# Allocate a powerful cloud machine on GCP
cloud_machine = inductiva.resources.MachineGroup(
    provider="GCP",
    machine_type="c3d-standard-360",
    spot=True
)

# Set up the OpenFOAM simulator
OpenFOAM = inductiva.simulators.OpenFOAM(
    version="2412",
    distribution="esi"
)

# Launch the simulation
task = OpenFOAM.run(
    input_dir="/path/to/highLiftConfiguration",
    shell_script="./Allrun",
    on=cloud_machine
)

# Monitor, download results and clean up
task.wait()
cloud_machine.terminate()
task.download_outputs()
task.print_summary()
```

>**Note**: Spot clusters are much cheaper than non spot, but the wholse simulation fails if one of the machines gets preempted.

## 📊 Simulation Summary

Once the simulation finishes, you’ll see output like this:

```bash
inductiva tasks info snkbbg86x8okt252gpryo535j

Task status: Success

Timeline:
	Waiting for Input         at 23/04, 11:21:02      2.696 s
	In Queue                  at 23/04, 11:21:05      58.487 s
	Preparing to Compute      at 23/04, 11:22:03      8.366 s
	In Progress               at 23/04, 11:22:12      147980.28 s
        ...
	Finalizing                at 25/04, 04:28:32      782.174 s
	Success                   at 25/04, 04:41:34      

Data:
	Size of zipped output:    25.00 GB
	Size of unzipped output:  31.55 GB
	Number of output files:   50712

Estimated computation cost (US$): 140.35 US$

Go to https://console.inductiva.ai/tasks/snkbbg86x8okt252gpryo535j for more details.
```

> The simulation ran for approximately **41 hours and 6 minutes**, as shown in the `In Progress` line.

## 🧩 What’s Next?

This benchmark shows how easily you can offload large-scale simulations to the
cloud. In the [next section](mpi-cluster-exafoam-high-lift-airfoil-benchmark), we’ll explore how to
**scale the simulation across multiple machines** to speed things up even further.

Stay tuned!
