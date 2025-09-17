# Run Your First Simulation
This tutorial will show you how to run FUNWAVE simulations using the Inductiva API. 

We will cover the `standing_waves` use case from the [BENCHMARK_FUNWAVE GitHub repository](https://github.com/fengyanshi/BENCHMARK_FUNWAVE/tree/master) to help you get started with simulations.


## Prerequisites

First, download the required files from the [repository](https://github.com/fengyanshi/BENCHMARK_FUNWAVE/tree/master/standing_waves).

After downloading, move all the files inside the `initial` folder into the `work` folder.

Your `work` folder should then look like this:

```
ls -lg
total 1632
-rw-rw-r--@ 1 staff     178 Nov 10  2024 Grid_Range.out
-rw-rw-r--@ 1 staff   29654 Nov 10  2024 LOG.txt
-rw-rw-r--@ 1 staff    3603 Nov 10  2024 eta0.txt
-rwxr-xr-x@ 1 staff  768584 Nov 10  2024 funwave_ab
-rw-rw-r--@ 1 staff    4987 Sep 17 11:11 input.txt
-rw-rw-r--@ 1 staff      17 Nov 10  2024 stat.txt
-rw-rw-r--@ 1 staff    1750 Nov 10  2024 time_dt.out
-rw-rw-r--@ 1 staff    3603 Nov 10  2024 u0.txt
-rw-rw-r--@ 1 staff    3603 Nov 10  2024 v0.txt
```

---

### Editing the `input.txt` file

Next, we need to update the configuration in `input.txt`:

1. **Set the output folder**
   Change the results directory to:

   ```
   RESULT_FOLDER = ./output/
   ```

2. **Fix the input file paths**
   Update the following variables to match the files you just copied into the `work` folder:

   ```
   ETA_FILE = eta0.txt
   U_FILE   = u0.txt
   V_FILE   = v0.txt
   ```

3. **Adjust CPU configuration**
   To run the simulation on 16 vCPUs, set:

   ```
   PX = 4
   PY = 4
   ```

That’s it! 🎉
Your setup is now ready, and you have everything needed to run the simulation.


## Running a FUNWAVE Simulation
Here is the code required to run a FUNWAVE simulation using the Inductiva API:

```python
import inductiva

machine_group = inductiva.resources.MachineGroup(
    provider="GCP",
    machine_type="c3d-highcpu-16",
    spot=True)

funwave = inductiva.simulators.FUNWAVE( \
    version="3.6")

task = funwave.run(
    input_dir="/Path/to/work/folder",
    sim_config_filename="input.txt",
    n_vcpus=16,
    on=machine_group)

# Wait for the simulation to finish and download the results
task.wait()
cloud_machine.terminate()

task.download_outputs()

task.print_summary()
```

> **Note**: Setting `spot=True` enables the use of [spot machines](../how-it-works/machines/spot-machines.md), which are available at substantial discounts. 
> However, your simulation may be interrupted if the cloud provider reclaims the machine.

To adapt this script for other FUNWAVE simulations, replace `input_dir` with the
path to your FUNWAVE input files and set the `case_name` accordingly.

When the simulation is complete, we terminate the machine, download the results and print a summary of the simulation as shown below.

```
Task status: Success

Timeline:
	Waiting for Input         at 17/09, 11:12:02      0.911 s
	In Queue                  at 17/09, 11:12:03      35.114 s
	Preparing to Compute      at 17/09, 11:12:38      9.047 s
	In Progress               at 17/09, 11:12:47      36.94 s
		├> 1.082 s         cp /FUNWAVE-TVD-Version_3.6/Makefile .
		├> 10.089 s        make
		├> 23.11 s         /opt/openmpi/4.1.6/bin/mpirun --np 16 --use-hwthread-cpus funwave-work/compiled_funwave input.txt
		├> 1.098 s         rm -r funwave-work
		└> 1.087 s         rm Makefile
	Finalizing                at 17/09, 11:13:24      0.59 s
	Success                   at 17/09, 11:13:24      

Data:
	Size of zipped output:    847.00 KB
	Size of unzipped output:  1.77 MB
	Number of output files:   216

Estimated computation cost (US$): 0.0021 US$
```

As you can see in the "In Progress" line, the part of the timeline that represents the actual execution of the simulation, 
the core computation time of this simulation was approximately 37 seconds.

```{banner_small}
:origin: funwave-quick-start
```

